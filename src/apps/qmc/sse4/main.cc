#include <alps/osiris/comm.h>

#include "sse.h"

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
	try {
#endif
	return alps::scheduler::start(argc, argv, 
							alps::scheduler::SimpleMCFactory<SSE_run>());
#ifndef BOOST_NO_EXCEPTIONS
	} catch (std::exception& e) {
	    std::cerr << e.what() << "\n";
		alps::comm_exit(true);
		return -1;
	} catch (...) {
		std::cerr << "Fatal Error: Unknown Exception!\n";
		return -2;
	}
#endif
}
