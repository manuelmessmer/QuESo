#pragma once

#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <string>

namespace queso {
namespace Testing {

	class GlobalConfig {
    public:
        static GlobalConfig& GetInstance() {
            static GlobalConfig config;
            return config;

        }
        std::string BaseDir;

    private:
        GlobalConfig() {
            auto& master = boost::unit_test::framework::master_test_suite();
            if (master.argc <= 1) {
                throw std::runtime_error(
                    "Test input directory not provided! "
                    "Usage: " + std::string(master.argv[0]) + " <input_dir>"
                );
            }
            BaseDir = master.argv[1];
        }
    };
} // End namespace Testing
} // End namespace queso
