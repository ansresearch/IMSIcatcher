//
// Copyright 2010-2012,2014 Ettus Research LLC
// Copyright 2018 Ettus Research, a National Instruments Company
//
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <uhd/exception.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/static.hpp>
#include <uhd/utils/thread.hpp>
#include <stdint.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/program_options.hpp>
#include <chrono>
#include <csignal>
#include <iostream>
#include <fstream>
#include <string>
#include <thread>

#include <stdlib.h>
#include <fftw3.h>


namespace po = boost::program_options;

/***********************************************************************
 * Signal handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

// clock counter
static volatile int tick = 0;

void clocker_func(double hoptime)
{
        while(1) {
                std::this_thread::sleep_for(
                        std::chrono::milliseconds(int64_t(hoptime))
                        );
                tick ++;
        }
}


// from https://home.zhaw.ch/kunr/NTM1/literatur/LTE%20in%20a%20Nutshell%20-%20Physical%20Layer.pdf
// What follows is for framing called "Normal CP":
// LTE subframe is transmitted in 0.5ms and is composed by seven OFDM symbols, each one made by
// N * 128 samples and a CP made of N * 10 (the first symbol) or N * 9 (the following six).
// A subframe is hence composed of N * (10 + 128 + 6 * (9 + 128)) = N * 960.
// By dividing number of samples by subframe duration we get Sampling rate is N * 1.92MS/s.
// When building the subframe, before IFFT only (N * 75 + 1) carriers are occupied
// and remaining (N * 53 - 1) are guards.
// N can be [1 2 4 8 12 16].
// For jamming efficiently we keep the same frame structure and we populate a contiguous
// subset of carriers randomly placed in the spectrum with QPSK symbols. We then change the
// position periodically.

// note: after IFFT the value of the signal in the temporal domain is bounded by the
// number of occupied carriers, i.e., if NT = 90, then assuming they are all
// populated with 1 + j, the signal is a sinc with amplitude limited to 90 * sqrt(2).
// So, if we plan to use short (sc16), we need to convert the fft_complex into short
// and multiply it by ~28000 / (90 * sqrt(2))
void create_lte_subframe (std::vector<std::complex<short> > &frame, int fft_length,
                          int subcarrier0_to_jam, int subcarriers_to_jam,
                          fftw_plan *p, fftw_complex *in, fftw_complex *out)
{
    // check it's valid fft length
    switch (fft_length) {
    case 128:
    case 256:
    case 512:
    case 1024:
    case 1536:
    case 2048:
        break;
    default:
        std::cerr << "Invalid fft size" << std::endl;
        exit (1);
    }

    // check number of elements in frame matches fft length
    int N = fft_length >> 7;
    int frame_length = N * (10 + 128 + 6 * (9 + 128));
    if (frame.size () != frame_length) {
        std::cerr << "frame is too short" << std::endl;
        exit (1);
    }

    // start filling subframe from head
    int current_symbol_position = 0;
    int ofdm_symbol_number = 0;

#define OFMD_SYMBOLS 7
    for (int ofdm_symbol_number = 0; ofdm_symbol_number < OFMD_SYMBOLS; ofdm_symbol_number ++) {
        // clean signal in frequency domain
        memset ((void *) in, 0, sizeof(fftw_complex) * fft_length);

        // now go around all carriers and fill qpsk coefficients
        for (int carrierjj = 0; carrierjj < subcarriers_to_jam; carrierjj ++) {
            int index = subcarrier0_to_jam - subcarriers_to_jam / 2 + carrierjj;
            index = (index + fft_length) % fft_length;
            in [index][0] = -1;
            in [index][1] = -1;
            if (drand48 () > 0.5) in [index][0] = 1;
            if (drand48 () > 0.5) in [index][1] = 1;
        }

        // get signal in time domain
        fftw_execute (*p);

        // copy and normalize ofdp symbol inside frame
        double normalization_factor = 28000 / (subcarriers_to_jam * sqrt(2));
        // fill CP first
        int cp_length = N * 9;
        if (ofdm_symbol_number == 0) cp_length = N * 10;
        int start_of_cp_in_symbol = fft_length - cp_length;
        for (int kk = 0; kk < cp_length; kk ++) {
            double real = out [start_of_cp_in_symbol + kk][0] * normalization_factor;
            double imag = out [start_of_cp_in_symbol + kk][1] * normalization_factor;
            std::complex<short> tmp (real, imag);
            frame [current_symbol_position + kk] = tmp;
        }
        for (int kk = 0; kk < fft_length; kk ++) {
            double real = out [kk][0] * normalization_factor;
            double imag = out [kk][1] * normalization_factor;
            std::complex<short> tmp (real, imag);
            frame [current_symbol_position + cp_length + kk] = tmp;
        }
        current_symbol_position = current_symbol_position + cp_length + fft_length;
    }
}

/***********************************************************************
 * Main function
 **********************************************************************/
int UHD_SAFE_MAIN(int argc, char* argv[])
{
    uhd::set_thread_priority_safe();

    // variables to be set by po
    std::string args, ant, subdev, ref, pps, otw, channel_list, freqs_list;
    size_t spb;
    double rate = 0;
    double gain, bw, lo_offset, hoptime;
    std::string file1, file2;
    std::string lte;
    int lte_frames, lte_N;
    int fft_length, lte_frame_length, subcarriers_to_jam;
    fftw_complex *in;
    fftw_complex *out;
    fftw_plan p;
    std::vector<std::complex<short> > lte_frame1;
    std::vector<std::complex<short> > lte_frame2;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("args", po::value<std::string>(&args)->default_value(""), "single uhd device address args")
        ("rate", po::value<double>(&rate), "rate of outgoing samples")
	("freqs", po::value<std::string>(&freqs_list)->default_value(""), "Frequencies to jam in Hz")
	("hoptime", po::value<double>(&hoptime)->default_value(1000.0), "hopping time in ms")
        ("lo-offset", po::value<double>(&lo_offset)->default_value(0.0),
            "Offset for frontend LO in Hz (optional)")
        ("gain", po::value<double>(&gain), "gain for the RF chain")
        ("ant", po::value<std::string>(&ant), "antenna selection")
        ("subdev", po::value<std::string>(&subdev), "subdevice specification")
        ("bw", po::value<double>(&bw), "analog frontend filter bandwidth in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "clock reference (internal, external, mimo, gpsdo)")
        ("pps", po::value<std::string>(&pps), "PPS source (internal, external, mimo, gpsdo)")
        ("otw", po::value<std::string>(&otw)->default_value("sc16"), "specify the over-the-wire sample mode")
        ("channels", po::value<std::string>(&channel_list)->default_value("0"), "which channels to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("int-n", "tune USRP with integer-N tuning")
        ("file1", po::value<std::string>(&file1), "file to transmit from channel 1")
        ("file2", po::value<std::string>(&file2), "file to transmit from channel 2")
	("lte", po::value<std::string>(&lte), "lte parameters")

    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("UHD TX Waveforms %s") % desc << std::endl;
        return ~0;
    }

    int files = 0;
    if (vm.count ("file1")) files ++;
    if (vm.count ("file2")) files ++;
    if (files == 1) {
	    throw std::runtime_error("If settings files, two must be set (file1 and file2)");
    }

    if (vm.count ("lte")) {
        if (files > 0) {
            throw std::runtime_error("lte settings not compatible with file1 and file2");
        }
        if (vm.count ("rate")) {
            throw std::runtime_error("rate setting not compatible with lte");
        }
        // for lte we expect separated by commas
        // - the N multiplier (1 2 4 8 12 16)
        // - the number of occupied carriers
        // - how many frames before changing central jamming carrier
        std::vector<std::string> lte_strings;
        boost::split(lte_strings, lte, boost::is_any_of("\"',"));
        if (lte_strings.size() != 3) {
            throw std::runtime_error("lte configuration is wrong: expect N,carrier2jam,frames");
        }
    
        lte_N = std::stoi(lte_strings[0]);
        switch (lte_N) {
            case 1:
            case 2:
            case 4:
            case 8:
            case 12:
            case 16:
                break;
            default:
                throw std::runtime_error("lte_N multiplier must be 1, 2, 4, 8, 12 or 16");
        }
        fft_length = lte_N * 128;
        lte_frame_length = lte_N * (10 + 128 + 6 * (9 + 128));

        subcarriers_to_jam = std::stoi(lte_strings[1]);
        if (subcarriers_to_jam < 1 || subcarriers_to_jam > fft_length) {
            throw std::runtime_error("number of carrier to jam is not valid");
        }

        lte_frames = std::stoi(lte_strings[2]);
        if (lte_frames < 1) {
            throw std::runtime_error("number of lte frames after changing frequency must be positive");
        }

        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_length);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_length);
        if (!in || !out) {
            throw std::runtime_error("cannot allocate memory for generating lte compliant frame");
        }
        p = fftw_plan_dft_1d(fft_length, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        lte_frame1.resize (lte_frame_length);
        lte_frame2.resize (lte_frame_length);

        rate = 1920000 * lte_N;
        std::cout << "lte configuration: fftsize = " << fft_length << ", frame_length = " << lte_frame_length << ", rate = " << rate << std::endl;
        std::cout << "lte configuration: number_of_carrier_to_jam = " << subcarriers_to_jam << std::endl;
        std::cout << "lte configuration: change central frequency every " << lte_frames << " half subframes" << std::endl;
    }

    // create a usrp device
    std::cout << std::endl;
    std::cout << boost::format("Creating the usrp device with: %s...") % args
              << std::endl;
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);
 
    // always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("subdev"))
        usrp->set_tx_subdev_spec(subdev);

    // detect which channels to use
    std::vector<std::string> channel_strings;
    std::vector<size_t> channel_nums;
    boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));
    if (channel_strings.size() != 2) {
        throw std::runtime_error("Need to have exactly two channels.");
    }

    for (size_t ch = 0; ch < channel_strings.size(); ch++) {
        size_t chan = std::stoi(channel_strings[ch]);
        if (chan >= usrp->get_tx_num_channels())
            throw std::runtime_error("Invalid channel(s) specified.");
        else
            channel_nums.push_back(std::stoi(channel_strings[ch]));
    }

    // Lock mboard clocks
    usrp->set_clock_source(ref);

    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

    // set the sample rate
    if (rate == 0) {
        std::cerr << "Please specify the sample rate with --rate" << std::endl;
        return ~0;
    }

    // extract frequencies
    std::vector<std::string> freqs_strings;
    std::vector<double> freqs;
    boost::split(freqs_strings, freqs_list, boost::is_any_of("\"',"));
    if (freqs_strings.size() < 1) {
        throw std::runtime_error("At least one frequency is needed.");
    }

    for (size_t fjj = 0; fjj < freqs_strings.size(); fjj ++) {
        freqs.push_back(std::stod(freqs_strings[fjj]));
    }

    for (size_t ch = 0; ch < channel_nums.size(); ch++) {
       std::cout << boost::format("Setting TX Rate: %f Msps...") % (rate / 1e6) << std::endl;
       usrp->set_tx_rate(rate, channel_nums[ch]);
       std::cout << boost::format("Actual TX Rate: %f Msps...") % (usrp->get_tx_rate(channel_nums[ch]) / 1e6)
                 << std::endl
                 << std::endl;

	std::cout << boost::format("Setting TX Freq1: %f MHz...") % (freqs[0] / 1e6)
	          << std::endl;
	std::cout << boost::format("Setting TX LO Offset: %f MHz...") % (lo_offset / 1e6)
	          << std::endl;
        uhd::tune_request_t tune_request(freqs[0], lo_offset);

        if (vm.count("int-n"))
            tune_request.args = uhd::device_addr_t("mode_n=integer");
        usrp->set_tx_freq(tune_request, channel_nums[ch]);
        std::cout << boost::format("Actual TX Freq: %f MHz...")
                         % (usrp->get_tx_freq(channel_nums[ch]) / 1e6)
                  << std::endl
                  << std::endl;

        // set the rf gain
        if (vm.count("gain")) {
            std::cout << boost::format("Setting TX Gain: %f dB...") % gain << std::endl;
            usrp->set_tx_gain(gain, channel_nums[ch]);
            std::cout << boost::format("Actual TX Gain: %f dB...")
                             % usrp->get_tx_gain(channel_nums[ch])
                      << std::endl
                      << std::endl;
        }

        // set the analog frontend filter bandwidth
        if (vm.count("bw")) {
            std::cout << boost::format("Setting TX Bandwidth: %f MHz...") % bw
                      << std::endl;
            usrp->set_tx_bandwidth(bw, channel_nums[ch]);
            std::cout << boost::format("Actual TX Bandwidth: %f MHz...")
                             % usrp->get_tx_bandwidth(channel_nums[ch])
                      << std::endl
                      << std::endl;
        }

        // set the antenna
        if (vm.count("ant"))
            usrp->set_tx_antenna(ant, channel_nums[ch]);
    }

    std::this_thread::sleep_for(std::chrono::seconds(1)); // allow for some setup time

    // create a transmit streamer
    // linearly map channels (index0 = channel0, index1 = channel1, ...)
    uhd::stream_args_t stream_args("sc16", otw);
    stream_args.channels             = channel_nums;
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);

    // create thread for counting time
    boost::thread_group clocker;
    clocker.create_thread(boost::bind(&clocker_func, hoptime));

    // load or precompute data
#define DEFAULT_SAMP    10240
    std::vector<std::complex<short>> buff(DEFAULT_SAMP);

    std::vector<std::complex<short>> samples1;
    std::vector<std::complex<short>> samples2;

    if (vm.count ("lte")) {
    }
    else if (files == 2) {
        std::ifstream infile1(file1.c_str (), std::ifstream::binary);
        std::ifstream infile2(file2.c_str (), std::ifstream::binary);

        while (!infile1.eof ()) {
            infile1.read ((char*)&buff.front(), buff.size()*sizeof(std::complex<short>));
            int read_bytes = infile1.gcount ();
            int read_elements = read_bytes / sizeof(std::complex<short>);
            if (read_elements < DEFAULT_SAMP) {
                    buff.resize (read_elements);
            }
            samples1.insert (samples1.end (), buff.begin (), buff.end ());
        }
        infile1.close();

        while (!infile2.eof ()) {
            infile2.read ((char*)&buff.front(), buff.size()*sizeof(std::complex<short>));
            int read_bytes = infile2.gcount ();
            int read_elements = read_bytes / sizeof(std::complex<short>);
            if (read_elements < DEFAULT_SAMP) {
                    buff.resize (read_elements);
            }
            samples2.insert (samples2.end (), buff.begin (), buff.end ());
        }
        infile2.close();

        std::cout << "file1 has " << samples1.size () << " elements" << std::endl;
        std::cout << "file2 has " << samples2.size () << " elements" << std::endl;

        if (samples1.size () != samples2.size ()) {
            std::cerr << "Files have different sizes" << std::endl;
            return 1;
        }
    } else {
	std::cout << "generating random samples" << std::endl;
#define STD_NOISE_LENGTH 1000000
	samples1.resize (STD_NOISE_LENGTH);
	samples2.resize (STD_NOISE_LENGTH);
	if (samples1.size () != STD_NOISE_LENGTH ||
	    samples2.size () != STD_NOISE_LENGTH) {
	    throw std::runtime_error("Cannot resize noise vectors");
	}
	for (int kk = 0; kk < STD_NOISE_LENGTH; kk ++) {
	    std::complex<short> tmp1 ((short) (2.0 * (drand48 () - 0.5) * 12000.0), (short) (2.0 * (drand48 () - 0.5) * 12000.0));
	    samples1[kk] = tmp1;
	    std::complex<short> tmp2 ((short) (2.0 * (drand48 () - 0.5) * 12000.0), (short) (2.0 * (drand48 () - 0.5) * 12000.0));
	    samples2[kk] = tmp2;

	    // std::cout << kk << ": " << samples1[kk] << " " << samples2[kk] << std::endl;
	}
    }

    size_t num_tx_samps = samples1.size ();

    std::vector<std::complex<short> *> banks;
    banks.push_back(&samples1.front ());
    banks.push_back(&samples2.front ());


    std::cout << boost::format("Setting device timestamp to 0...") << std::endl;
    if (channel_nums.size() > 1) {
        // Sync times
        if (pps == "mimo") {
            UHD_ASSERT_THROW(usrp->get_num_mboards() == 2);

            // make mboard 1 a slave over the MIMO Cable
            usrp->set_time_source("mimo", 1);

            // set time on the master (mboard 0)
            usrp->set_time_now(uhd::time_spec_t(0.0), 0);

            // sleep a bit while the slave locks its time to the master
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        } else {
            if (pps == "internal" or pps == "external" or pps == "gpsdo")
                usrp->set_time_source(pps);
            usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
            std::this_thread::sleep_for(
                std::chrono::seconds(1)); // wait for pps sync pulse
        }
    } else {
        usrp->set_time_now(0.0);
    }

    // Check Ref and LO Lock detect
    std::vector<std::string> sensor_names;
    const size_t tx_sensor_chan = channel_nums.empty() ? 0 : channel_nums[0];
    sensor_names                = usrp->get_tx_sensor_names(tx_sensor_chan);
    if (std::find(sensor_names.begin(), sensor_names.end(), "lo_locked")
        != sensor_names.end()) {
        uhd::sensor_value_t lo_locked = usrp->get_tx_sensor("lo_locked", tx_sensor_chan);
        std::cout << boost::format("Checking TX: %s ...") % lo_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(lo_locked.to_bool());
    }
    const size_t mboard_sensor_idx = 0;
    sensor_names                   = usrp->get_mboard_sensor_names(mboard_sensor_idx);
    if ((ref == "mimo")
        and (std::find(sensor_names.begin(), sensor_names.end(), "mimo_locked")
                != sensor_names.end())) {
        uhd::sensor_value_t mimo_locked =
            usrp->get_mboard_sensor("mimo_locked", mboard_sensor_idx);
        std::cout << boost::format("Checking TX: %s ...") % mimo_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(mimo_locked.to_bool());
    }
    if ((ref == "external")
        and (std::find(sensor_names.begin(), sensor_names.end(), "ref_locked")
                != sensor_names.end())) {
        uhd::sensor_value_t ref_locked =
            usrp->get_mboard_sensor("ref_locked", mboard_sensor_idx);
        std::cout << boost::format("Checking TX: %s ...") % ref_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(ref_locked.to_bool());
    }

    std::signal(SIGINT, &sig_int_handler);
    std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

    // Set up metadata. We start streaming a bit in the future
    // to allow MIMO operation:
    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst   = false;
    md.has_time_spec  = true;
    md.time_spec      = usrp->get_time_now() + uhd::time_spec_t(0.1);

    // send data until the signal handler gets called
    // or if we accumulate the number of samples specified (unless it's 0)
    uint64_t num_acc_samps = 0;
    int tick_next = tick + 2;
    int freqjj = 1;

    int lte_count_frames = 0;
    int subcarrier0_to_jam_chain1 = 0;
    int subcarrier0_to_jam_chain2 = 0;
    int regenerated = 0;
    while (true) {
        // Break on the end of duration or CTRL-C
        if (stop_signal_called) {
            break;
        }

        // send the entire contents of the buffer
        if (vm.count ("lte")) {
            if (lte_count_frames == 0) {
                regenerated ++;
                subcarrier0_to_jam_chain1 = subcarriers_to_jam / 2 + (fft_length - subcarriers_to_jam + 1) * drand48 () - fft_length / 2;
                subcarrier0_to_jam_chain2 = subcarriers_to_jam / 2 + (fft_length - subcarriers_to_jam + 1) * drand48 () - fft_length / 2;
                create_lte_subframe (lte_frame1, fft_length,
                                     subcarrier0_to_jam_chain1, subcarriers_to_jam,
                                     &p, in, out);
                create_lte_subframe (lte_frame2, fft_length,
                                     subcarrier0_to_jam_chain2, subcarriers_to_jam,
                                     &p, in, out);
                lte_count_frames = lte_frames;
            }
            lte_count_frames --;

            num_tx_samps = lte_frame1.size ();
            std::vector<std::complex<short> *> ltebanks;
            ltebanks.push_back(&lte_frame1.front ());
            ltebanks.push_back(&lte_frame2.front ());
            num_acc_samps += tx_stream->send(ltebanks, num_tx_samps, md);
        } else {
            num_acc_samps += tx_stream->send(banks, num_tx_samps, md);
        }

        md.start_of_burst = false;
        md.has_time_spec  = false;

	if (freqs.size () > 1 && tick >= tick_next) {
            tick_next = tick + 1;

            double freq = freqs[freqjj];
            uhd::tune_request_t tune_request;
            std::cout << boost::format("Setting TX Freq: %f MHz...") % (freq/1e6) << std::endl;
            uhd::tune_request_t tune_request_tmp(freq, lo_offset);
            if (vm.count("int-n"))
                tune_request_tmp.args = uhd::device_addr_t("mode_n=integer");
            usrp->set_tx_freq(tune_request_tmp, 0);
            usrp->set_tx_freq(tune_request_tmp, 1);
            freqjj ++;
            if (freqjj == freqs.size ()) freqjj = 0;

            std::cout << boost::format("Actual TX Freq0: %f MHz...") % (usrp->get_tx_freq(0)/1e6);
            std::cout << " txrate = " << usrp->get_tx_rate(channel_nums[0]) / 1e6  << std::endl << std::endl;
            std::cout << boost::format("Actual TX Freq1: %f MHz...") % (usrp->get_tx_freq(1)/1e6);
            std::cout << " txrate = " << usrp->get_tx_rate(channel_nums[0]) / 1e6  << std::endl << std::endl;
            if (vm.count ("lte")) {
                std::cout << "lte frame regenerated " << regenerated << " times" << std::endl;
                regenerated = 0;
            }

	}
    }

    // send a mini EOB packet
    md.end_of_burst = true;
    tx_stream->send("", 0, md);

    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;
    return EXIT_SUCCESS;
}
