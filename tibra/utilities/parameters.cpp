#include "utilities/parameters.h"

namespace tibra {


std::ostream& operator<< (std::ostream& rOStream, const Parameters& rThis){
    rThis.PrintInfo(rOStream);
    return rOStream;
}

} // End tibra namespace