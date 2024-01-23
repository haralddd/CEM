//
// Created by haralg on 15.12.23.
//

#ifndef CEM_PARSE_ARGS_H
#define CEM_PARSE_ARGS_H

#include "args-parser/all.hpp"

int parse_args(int argc, char **argv) {
    using namespace Args;

    try {
        CmdLine cmd(argc, argv);

        Arg input_file('i', true, true);
        input_file.setDescription("Input file name");
        input_file.setLongDescription("Input file name");

        Arg output_file('o', true, true);
        output_file.setDescription("Output file name");
        output_file.setLongDescription("Output file name");

        Arg generate_file('g', false, false);
        generate_file.setDescription("Generate a surface");
        generate_file.setLongDescription("Generate a surface");

        Help help;
        help.setExecutable(argv[0]);
        help.setDescription("This program computes the correlation function of a surface");

        cmd.addArg(input_file);
        cmd.addArg(output_file);
        cmd.addArg(generate_file);
        cmd.addArg(help);
    } catch (const std::exception &e) {
        std::cerr << "Argument parsing failed with error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}


#endif //CEM_PARSE_ARGS_H
