#include "machineInfo.h"

#ifdef WIN32
#	include <windows.h>
#else
#	include <sys/sysctl.h>
#endif

#define ARRAY_SIZE(a) (sizeof (a) / sizeof ((a)[0]))


int StMachineInfo::numProcs(void) {

#	ifdef WIN32
	int numCPU;
	SYSTEM_INFO sysinfo;
	GetSystemInfo( &sysinfo );
	numCPU = sysinfo.dwNumberOfProcessors;
	return numCPU;

#	else

	int numCPU = 0;
	int nprocs;
	size_t len = sizeof(nprocs); 
	static int mib[2] = { CTL_HW, HW_NCPU };

	/* get the number of CPUs from the system */
	sysctl(mib, 2, &numCPU, &len, NULL, 0);

	if( numCPU < 1 ) 
		{
		mib[1] = HW_NCPU;

		if (sysctl (mib, ARRAY_SIZE(mib), &nprocs, &len, NULL, 0) == 0 && len == sizeof (nprocs) && 0 < nprocs)
			numCPU = nprocs;

		if( numCPU < 1 )
			numCPU = 1;
		}
	return numCPU;
#	endif
}










