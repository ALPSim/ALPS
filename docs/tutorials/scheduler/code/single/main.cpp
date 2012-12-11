// Copyright (C) 20012 TODO: Licence needet!

#include <ising.hpp>

int main(int argc, char *argv[]) {

    try {
        args options(argc, argv);

        alps::parameters_type<ising_sim>::type parameters;

        std::string suffix = options.inputfile.substr(options.inputfile.find_last_of('.'));
        if (suffix == ".xml")
            parameters = alps::make_parameters_from_xml(options.inputfile);
        else if (suffix == ".h5")
            alps::hdf5::archive(options.inputfile)["/parameters"] >> parameters;
        else
            throw std::runtime_error("Unsupported input format: " + suffix + "!");

        ising_sim sim(parameters);

        if (options.resume && boost::filesystem::exists(options.checkpointfile))
            sim.load(options.checkpointfile);

        sim.run(stop_callback(options.timelimit));

                // make checkpoint
        sim.save(options.checkpointfile);

        using alps::collect_results;
        alps::results_type<ising_sim>::type results = collect_results(sim);

        std::cout << results << std::endl;
        alps::hdf5::archive ar(options.outputfile, "w");
        ar["/parameters"] << parameters;
        ar["/simulation/results"] << results;

    } catch (std::exception const & e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Caught unknown exception" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
