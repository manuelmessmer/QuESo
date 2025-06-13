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

#include <string>
#include <filesystem>

namespace queso {
namespace Testing {

///@name QuESo classes
///@{

/// @class  TemporaryFile.
/// @author Manuel Messmer
/// @brief  Allows to create temporary files. Removes files from disc on destruction.
class TemporaryFile {
public:
    ///@name Life cycle
    ///@{

    /// Constructor
    explicit TemporaryFile(const std::filesystem::path& rFilename)
        : mFilename(std::filesystem::temp_directory_path() / rFilename)
    {}

    /// Destructor
    ~TemporaryFile() {
        std::error_code ec;
        std::filesystem::remove(mFilename, ec);
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns filename.
    /// @return std::string
    std::string GetString() const {
        return mFilename.string();
    }

    /// @brief Returns filename.
    /// @return const std::filesystem::path&
    const std::filesystem::path& GetPath() const {
        return mFilename;
    }

private:
    ///@}
    ///@name Private members
    ///@{

    std::filesystem::path mFilename;
    ///@}
};
///@}

} // End Testing namespace
} // End queso namespace