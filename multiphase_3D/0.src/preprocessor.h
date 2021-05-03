/*use gpu*/
!#define _openacc
/*enable profiling gpu code*/
!#define gpu_profiling

/*mpi communication directions*/
/* x direction domain decomposition disabled */
#define y_mpi
#define z_mpi

/*old mrt_bb_opt=1, mrt_ori=2, srt=3, adv_opt=4*/
#define mrt 2

