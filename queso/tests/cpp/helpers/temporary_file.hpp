//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#include <atomic>
#include <chrono>
#include <filesystem>
#include <string>

namespace queso::testing {

    ///@name QuESo classes
    ///@{

    /// @class  TemporaryFile.
    /// @author Manuel Messmer
    /// @brief  Allows to create temporary files. Removes files from disc on destruction.
    class TemporaryFile
    {
    public:
        ///@name Life cycle
        ///@{

        /// Constructor
        explicit TemporaryFile(const std::filesystem::path& rFilename)
        {
            static std::atomic<std::size_t> sequence{ 0 };
            const auto timestamp = std::chrono::steady_clock::now().time_since_epoch().count();
            for (;;) {
                mDirectory =
                    std::filesystem::temp_directory_path()
                    / ("queso-test-" + std::to_string(timestamp) + "-" + std::to_string(sequence.fetch_add(1)));
                std::error_code error;
                if (std::filesystem::create_directory(mDirectory, error)) { break; }
                if (error) {
                    throw std::filesystem::filesystem_error("Cannot create test directory", mDirectory, error);
                }
            }
            mFilename = mDirectory / rFilename;
        }

        /// Destructor
        ~TemporaryFile()
        {
            std::error_code ec;
            std::filesystem::remove(mFilename, ec);
            std::filesystem::remove(mDirectory, ec);
        }

        TemporaryFile(const TemporaryFile&) = delete;
        TemporaryFile& operator=(const TemporaryFile&) = delete;

        ///@}
        ///@name Operations
        ///@{

        /// @brief Returns filename.
        /// @return std::string
        std::string GetString() const
        { return mFilename.string(); }

        /// @brief Returns filename.
        /// @return const std::filesystem::path&
        const std::filesystem::path& GetPath() const
        { return mFilename; }

    private:
        ///@}
        ///@name Private members
        ///@{

        std::filesystem::path mDirectory;
        std::filesystem::path mFilename;
        ///@}
    };
    ///@}

} // namespace queso::testing
