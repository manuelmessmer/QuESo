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

#ifndef BINARY_BUFFER_WRITER_HPP
#define BINARY_BUFFER_WRITER_HPP

//// STL includes
#include <fstream>
#include <vector>
#include <cstring>

//// Project includes
#include "queso/includes/define.hpp"

namespace queso {

///@name QuESo Classes
///@{

/// @class  BinaryBufferWriter
/// @author Manuel Messmer
/// @brief  Buffer to write data in binary format to file.
class BinaryBufferWriter {
public:
    ///@name Type definitions
    ///@{

    enum class EndianType {little, big};

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor.
    /// @param rOut Filestream.
    /// @param Endian Options: {big, little}.
    /// @param BufferSize (Default: 1MB).
    BinaryBufferWriter(std::ofstream& rOut, EndianType Endian, IndexType BufferSize = (int)1 << 20) : // default 1MB
          mOut(rOut),
          mEndian(Endian),
          mBufferSize(BufferSize),
          mBufferPos(0),
          mBuffer(BufferSize)
    {
    }

    /// Destructor
    ~BinaryBufferWriter() {
        Flush();
    }

    /// Delete copy and move operators.
    BinaryBufferWriter(const BinaryBufferWriter& rOther) = delete;
    BinaryBufferWriter& operator=(const BinaryBufferWriter& rOther) = delete;
    BinaryBufferWriter(BinaryBufferWriter&& rOther) = delete;
    BinaryBufferWriter& operator=(BinaryBufferWriter&& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Writes raw data.
    /// @param pData
    /// @param Size
    void WriteRaw(const char* pData, IndexType Size){
        AppendToBuffer(pData, Size);
    }

    /// @brief Writes basic values.
    /// @tparam TType
    /// @param rValue
    template<typename TType>
    void WriteValue(const TType& rValue) {
        static_assert(std::is_trivially_copyable_v<TType>, "T must be trivially copyable");
        TType temp = rValue;
        if(NativeEndian() != mEndian) {
            BinaryBufferWriter::SwapEnd(temp);
        }
        AppendToBuffer(reinterpret_cast<const char*>(&temp), sizeof(TType));
    }

    /// @brief Writes C-Style arrays.
    /// @tparam TType
    /// @tparam N
    /// @param rArray
    template<typename TType, IndexType N>
    void WriteArray(const TType (&rArray)[N]) {
        static_assert(std::is_trivially_copyable_v<TType>, "T must be trivially copyable");
        for (IndexType i = 0; i < N; ++i) {
            WriteValue(rArray[i]);
        }
    }

    /// @brief Flushes the buffer to the filestream.
    void Flush() {
        if (mBufferPos > 0) {
            mOut.write(mBuffer.data(), static_cast<long int>(mBufferPos));
            QuESo_ERROR_IF(!mOut) << "Write failed during Flush." << std::endl;
            mBufferPos = 0; // Reset buffer.
        }
    }

private:

    ///@}
    ///@name Private Operations
    ///@{

    /// @brief Adds given data to buffer.
    /// @param pData
    /// @param Size
    void AppendToBuffer(const char* pData, IndexType Size) {
        if(Size > mBufferSize) {
            Flush();
            // Write large data directly to stream (no buffering).
            mOut.write(pData, static_cast<long int>(Size));
            QuESo_ERROR_IF(!mOut) << "Write failed during Flush." << std::endl;
            return;
        }
        if (mBufferPos + Size > mBufferSize) {
            Flush();
        }
        std::memcpy(mBuffer.data() + mBufferPos, pData, Size);
        mBufferPos += Size;
    }

    /// @brief Returns the endian type of the current system.
    ///        This function should only be called once.
    ///        Please, use NativeEndian.
    /// @return EndianType.
    static EndianType DetectEndianOnce() {
        union {
            uint32_t i;
            char c[4];
        } test = {0x01020304};

        return (test.c[0] == 0x04) ? EndianType::little : EndianType::big;
    }

    /// @brief Returns the native endian type of the given system.
    /// @return EndianType.
    static EndianType NativeEndian() {
        // Cache value
        static const EndianType endian = DetectEndianOnce();
        return endian;
    }

    /// @brief Swaps bytes.
    /// @tparam TType
    /// @param rValue
    template<typename TType>
    static void SwapEnd(TType& rValue) {
        auto* bytes = reinterpret_cast<std::uint8_t*>(&rValue);
        IndexType n = sizeof(TType);
        for (IndexType i = 0; i < n / 2; ++i) {
            std::swap(bytes[i], bytes[n - 1 - i]);
        }
    }

    ///@}
    ///@name Member variables
    ///@{

    std::ofstream& mOut;
    EndianType mEndian;
    IndexType mBufferSize;
    IndexType mBufferPos;
    std::vector<char> mBuffer;

  ///@}
}; // End class BinaryBufferWriter
///@} End QuESo Classes

} // End namespace queso

#endif // BINARY_BUFFER_WRITER_HPP
