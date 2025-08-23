
// If 1, the scenes are not deleted after the rendering. Makes the RTPP faster, but consumes more memory
#define KEEP_SCENES_IN_MEMORY 0

// If 1, the dataset providers are not deleted after the rendering. Makes the RTPP faster, but consumes more memory
#define KEEP_PROVIDERS_IN_MEMORY 0

// If 1, the time to process and render is measured and the results are printed in the end
#define BASIC_PROFILE 1

// Choose the compositer to use for the transparencies
#define ICE_T_COMPOSITER 1
#define ALL_GATHER_TRANSPARENT_COMPOSITER 2
#define PARALLEL_DEPTH_PEEL_COMPOSITER 3

#define TRANSPARENT_COMPOSITER_TO_USE PARALLEL_DEPTH_PEEL_COMPOSITER
