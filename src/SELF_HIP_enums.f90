MODULE SELF_HIP_enums

    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ==============================================================================
    ! hipfort: FORTRAN INTERFACEs for GPU kernels
    ! ==============================================================================
    ! Copyright (c) 2020-2022 Advanced Micro Devices, Inc. All rights reserved.
    ! [MITx11 License]
    !
    ! Permission is hereby granted, free of charge, to any person obtaining a copy
    ! of this software and associated documentation files (the "Software"), to deal
    ! in the Software without restriction, including without limitation the rights
    ! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    ! copies of the Software, and to permit persons to whom the Software is
    ! furnished to do so, subject to the following conditions:
    !
    ! The above copyright notice and this permission notice shall be included in
    ! all copies or substantial portions of the Software.
    !
    ! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    ! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    ! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
    ! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    ! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    ! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    ! THE SOFTWARE.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
      USE ISO_C_BINDING
    
      IMPLICIT NONE
    
      !> Derived TYPE that can be mapped directly to a CUDA/HIP C++ dim3.
      TYPE,BIND(c) :: dim3
        INTEGER(C_INT) :: x = 1,y = 1,z = 1
      END TYPE dim3
    
      TYPE,BIND(c) :: hipDeviceProp_t ! as of ROCm 4.4
        CHARACTER(kind=C_CHAR) :: name(256)            !< Device name.
        INTEGER(C_SIZE_T) :: totalGlobalMem     !< Size of global memory region (in bytes).
        INTEGER(C_SIZE_T) :: sharedMemPerBlock  !< Size of shared memory region (in bytes).
        INTEGER(C_INT) :: regsPerBlock          !< Registers per block.
        INTEGER(C_INT) :: warpSize              !< Warp size.
        INTEGER(C_INT) :: maxThreadsPerBlock    !< Max work items per work group or workgroup max size.
        INTEGER(C_INT) :: maxThreadsDim(3)      !< Max number of threads in each dimension (XYZ) of a block.
        INTEGER(C_INT) :: maxGridSize(3)        !< Max grid dimensions (XYZ).
        INTEGER(C_INT) :: clockRate             !< Max clock frequency of the multiProcessors in khz.
        INTEGER(C_INT) :: memoryClockRate       !< Max global memory clock frequency in khz.
        INTEGER(C_INT) :: memoryBusWidth        !< Global memory bus width in bits.
        INTEGER(C_SIZE_T) :: totalConstMem      !< Size of shared memory region (in bytes).
        INTEGER(C_INT) :: major  !< Major compute capability.  On HCC, this is an approximation and features may
        !< differ from CUDA CC.  See the arch feature flags for portable ways to query
        !< feature caps.
        INTEGER(C_INT) :: minor  !< Minor compute capability.  On HCC, this is an approximation and features may
        !< differ from CUDA CC.  See the arch feature flags for portable ways to query
        !< feature caps.
        INTEGER(C_INT) :: multiProcessorCount          !< Number of multi-processors (compute units).
        INTEGER(C_INT) :: l2CacheSize                  !< L2 cache size.
        INTEGER(C_INT) :: maxThreadsPerMultiProcessor  !< Maximum resident threads per multi-processor.
        INTEGER(C_INT) :: computeMode                  !< Compute mode.
        INTEGER(C_INT) :: clockInstructionRate  !< Frequency in khz of the timer used by the device-side "clock*"
        !< instructions.  New for HIP.
        INTEGER(C_INT) arch       !< Architectural feature flags.  New for HIP.
        INTEGER(C_INT) :: concurrentKernels     !< Device can possibly execute multiple kernels concurrently.
        INTEGER(C_INT) :: pciDomainID           !< PCI Domain ID
        INTEGER(C_INT) :: pciBusID              !< PCI Bus ID.
        INTEGER(C_INT) :: pciDeviceID           !< PCI Device ID.
        INTEGER(C_SIZE_T) :: maxSharedMemoryPerMultiProcessor  !< Maximum Shared Memory Per Multiprocessor.
        INTEGER(C_INT) :: isMultiGpuBoard                      !< 1 if device is on a multi-GPU board, 0 if not.
        INTEGER(C_INT) :: canMapHostMemory                     !< Check whether HIP can map host memory
        INTEGER(C_INT) :: gcnArch                              !< DEPRECATED: use gcnArchName instead
        CHARACTER(kind=C_CHAR) :: gcnArchName(256)                    !< AMD GCN Arch Name.
        INTEGER(C_INT) :: integrated            !< APU vs dGPU
        INTEGER(C_INT) :: cooperativeLaunch            !< HIP device supports cooperative launch
        INTEGER(C_INT) :: cooperativeMultiDeviceLaunch !< HIP device supports cooperative launch on multiple devices
        INTEGER(C_INT) :: maxTexture1DLinear    !< Maximum size for 1D textures bound to linear memory
        INTEGER(C_INT) :: maxTexture1D          !< Maximum number of elements in 1D images
        INTEGER(C_INT) :: maxTexture2D(2)       !< Maximum dimensions (width, height) of 2D images, in image elements
        INTEGER(C_INT) :: maxTexture3D(3)       !< Maximum dimensions (width, height, depth) of 3D images, in image elements
        TYPE(C_PTR) :: hdpMemFlushCntl      !< Addres of HDP_MEM_COHERENCY_FLUSH_CNTL register
        TYPE(C_PTR) :: hdpRegFlushCntl      !< Addres of HDP_REG_COHERENCY_FLUSH_CNTL register
        INTEGER(C_SIZE_T) :: memPitch                 !<Maximum pitch in bytes allowed by memory copies
        INTEGER(C_SIZE_T) :: textureAlignment         !<Alignment requirement for textures
        INTEGER(C_SIZE_T) :: texturePitchAlignment    !<Pitch alignment requirement for texture references bound to pitched memory
        INTEGER(C_INT) :: kernelExecTimeoutEnabled    !<Run time limit for kernels executed on the device
        INTEGER(C_INT) :: ECCEnabled                  !<Device has ECC support enabled
        INTEGER(C_INT) :: tccDriver                   !< 1:If device is Tesla device using TCC driver, else 0
        INTEGER(C_INT) :: cooperativeMultiDeviceUnmatchedFunc        !< HIP device supports cooperative launch on multiple
        !devices with unmatched FUNCTIONs
        INTEGER(C_INT) :: cooperativeMultiDeviceUnmatchedGridDim     !< HIP device supports cooperative launch on multiple
        !devices with unmatched grid dimensions
        INTEGER(C_INT) :: cooperativeMultiDeviceUnmatchedBlockDim    !< HIP device supports cooperative launch on multiple
        !devices with unmatched block dimensions
        INTEGER(C_INT) :: cooperativeMultiDeviceUnmatchedSharedMem   !< HIP device supports cooperative launch on multiple
        !devices with unmatched shared memories
        INTEGER(C_INT) :: isLargeBar                  !< 1: if it is a large PCI bar device, else 0
        INTEGER(C_INT) :: asicRevision                !< Revision of the GPU in this device
        INTEGER(C_INT) :: managedMemory               !< Device supports allocating managed memory on this system
        INTEGER(C_INT) :: directManagedMemAccessFromHost !< Host can directly access managed memory on the device without migration
        INTEGER(C_INT) :: concurrentManagedAccess     !< Device can coherently access managed memory concurrently with the CPU
        INTEGER(C_INT) :: pageableMemoryAccess        !< Device supports coherently accessing pageable memory
        !< without calling hipHostRegister on it
        INTEGER(C_INT) :: pageableMemoryAccessUsesHostPageTables !< Device accesses pageable memory via the host's page tables
        CHARACTER(kind=C_CHAR) :: GPUFORT_PADDING(256) !< GPUFORT :Some extra bytes to prevent seg faults in newer versions
      END TYPE hipDeviceProp_t
    
      ! runtime api parameters
    
      INTEGER,PARAMETER :: hipIpcMemLazyEnablePeerAccess = 0
    
      INTEGER,PARAMETER :: hipStreamDefault = &
                           0  !< Default stream creation flags. These are used with hipStreamCreate = ().
      INTEGER,PARAMETER :: hipStreamNonBlocking = 1  !< Stream does not implicitly synchronize with null stream
    
      INTEGER,PARAMETER :: hipEventDefault = 0  !< Default flags
      INTEGER,PARAMETER :: hipEventBlockingSync = &
                           1  !< Waiting will yield CPU.  Power-friENDly and usage-friENDly but may increase latency.
      INTEGER,PARAMETER :: hipEventDisableTiming = &
                           2  !< Disable event's capability to record timing information.  May improve performance.
      INTEGER,PARAMETER :: hipEventInterprocess = 4  !< Event can support IPC.  @warning - not supported in HIP.
      INTEGER,PARAMETER :: hipEventReleaseToDevice = &
                           1073741824 !< 0x40000000 - Use a device-scope release when recording this event.  This flag is useful to
      !INTEGER, parameter :: hipEventReleaseToSystem = &
      !    2147483648 !< 0x80000000 - Use a system-scope release that when recording this event.  This flag is
      INTEGER,PARAMETER :: hipHostMallocDefault = 0
      INTEGER,PARAMETER :: hipHostMallocPortable = 1  !< Memory is considered allocated by all contexts.
      INTEGER,PARAMETER :: hipHostMallocMapped = &
                           2  !< Map the allocation into the address space for the current device.  The device pointer
      INTEGER,PARAMETER :: hipHostMallocWriteCombined = 4
      INTEGER,PARAMETER :: hipHostMallocNumaUser = &
                           536870912 !< 0x20000000 - Host memory allocation will follow numa policy set by user
      INTEGER,PARAMETER :: hipHostMallocCoherent = &
                           1073741824 !< 0x40000000 - Allocate coherent memory. Overrides HIP_COHERENT_HOST_ALLOC for specific
      !INTEGER, parameter :: hipHostMallocNonCoherent = &
      !    2147483648 !< 0x80000000 - Allocate non-coherent memory. Overrides HIP_COHERENT_HOST_ALLOC for specific
      INTEGER,PARAMETER :: hipMemAttachGlobal = 1    !< Memory can be accessed by any stream on any device
      INTEGER,PARAMETER :: hipMemAttachHost = 2    !< Memory cannot be accessed by any stream on any device
      INTEGER,PARAMETER :: hipMemAttachSingle = 4    !< Memory can only be accessed by a single stream on
      !< the associated device
      INTEGER,PARAMETER :: hipDeviceMallocDefault = 0
      INTEGER,PARAMETER :: hipDeviceMallocFinegrained = 1  !< Memory is allocated in fine grained region of device.
    
      INTEGER,PARAMETER :: hipHostRegisterDefault = 0   !< Memory is Mapped and Portable
      INTEGER,PARAMETER :: hipHostRegisterPortable = 1  !< Memory is considered registered by all contexts.
      INTEGER,PARAMETER :: hipHostRegisterMapped = &
                           2  !< Map the allocation into the address space for the current device.  The device pointer
      INTEGER,PARAMETER :: hipHostRegisterIoMemory = 4  !< Not supported.
      INTEGER,PARAMETER :: hipExtHostRegisterCoarseGrained = 8  !< Coarse Grained host memory lock
    
      INTEGER,PARAMETER :: hipDeviceScheduleAuto = 0  !< Automatically select between Spin and Yield
      INTEGER,PARAMETER :: hipDeviceScheduleSpin = &
                           1  !< Dedicate a CPU core to spin-wait.  Provides lowest latency, but burns a CPU core and
      INTEGER,PARAMETER :: hipDeviceScheduleYield = &
                           2  !< Yield the CPU to the operating system when waiting.  May increase latency, but lowers
      INTEGER,PARAMETER :: hipDeviceScheduleBlockingSync = 4
      INTEGER,PARAMETER :: hipDeviceScheduleMask = 7
    
      INTEGER,PARAMETER :: hipDeviceMapHost = 8
      INTEGER,PARAMETER :: hipDeviceLmemResizeToMax = 22 ! 16
    
      INTEGER,PARAMETER :: hipArrayDefault = 0  !< Default HIP array allocation flag
      INTEGER,PARAMETER :: hipArrayLayered = 1
      INTEGER,PARAMETER :: hipArraySurfaceLoadStore = 2
      INTEGER,PARAMETER :: hipArrayCubemap = 4
      INTEGER,PARAMETER :: hipArrayTextureGather = 8
    
      INTEGER,PARAMETER :: hipOccupancyDefault = 0
    
      INTEGER,PARAMETER :: hipCooperativeLaunchMultiDeviceNoPreSync = 1
      INTEGER,PARAMETER :: hipCooperativeLaunchMultiDeviceNoPostSync = 2
    
      ENUM,BIND(c)
        ENUMERATOR :: HIP_SUCCESS = 0
        ENUMERATOR :: HIP_ERROR_INVALID_VALUE
        ENUMERATOR :: HIP_ERROR_NOT_INITIALIZED
        ENUMERATOR :: HIP_ERROR_LAUNCH_OUT_OF_RESOURCES
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipMemoryTYPEHost
        ENUMERATOR :: hipMemoryTYPEDevice
        ENUMERATOR :: hipMemoryTYPEArray
        ENUMERATOR :: hipMemoryTYPEUnified
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipSuccess = 0
        ENUMERATOR :: hipErrorInvalidValue = 1
        ENUMERATOR :: hipErrorOutOfMemory = 2
        ENUMERATOR :: hipErrorMemoryAllocation = 2
        ENUMERATOR :: hipErrorNotInitialized = 3
        ENUMERATOR :: hipErrorInitializationError = 3
        ENUMERATOR :: hipErrorDeinitialized = 4
        ENUMERATOR :: hipErrorProfilerDisabled = 5
        ENUMERATOR :: hipErrorProfilerNotInitialized = 6
        ENUMERATOR :: hipErrorProfilerAlreadyStarted = 7
        ENUMERATOR :: hipErrorProfilerAlreadyStopped = 8
        ENUMERATOR :: hipErrorInvalidConfiguration = 9
        ENUMERATOR :: hipErrorInvalidPitchValue = 12
        ENUMERATOR :: hipErrorInvalidSymbol = 13
        ENUMERATOR :: hipErrorInvalidDevicePointer = 17
        ENUMERATOR :: hipErrorInvalidMemcpyDirection = 21
        ENUMERATOR :: hipErrorInsufficientDriver = 35
        ENUMERATOR :: hipErrorMissingConfiguration = 52
        ENUMERATOR :: hipErrorPriorLaunchFailure = 53
        ENUMERATOR :: hipErrorInvalidDeviceFUNCTION = 98
        ENUMERATOR :: hipErrorNoDevice = 100
        ENUMERATOR :: hipErrorInvalidDevice = 101
        ENUMERATOR :: hipErrorInvalidImage = 200
        ENUMERATOR :: hipErrorInvalidContext = 201
        ENUMERATOR :: hipErrorContextAlreadyCurrent = 202
        ENUMERATOR :: hipErrorMapFailed = 205
        ENUMERATOR :: hipErrorMapBufferObjectFailed = 205
        ENUMERATOR :: hipErrorUnmapFailed = 206
        ENUMERATOR :: hipErrorArrayIsMapped = 207
        ENUMERATOR :: hipErrorAlreadyMapped = 208
        ENUMERATOR :: hipErrorNoBinaryForGpu = 209
        ENUMERATOR :: hipErrorAlreadyAcquired = 210
        ENUMERATOR :: hipErrorNotMapped = 211
        ENUMERATOR :: hipErrorNotMappedAsArray = 212
        ENUMERATOR :: hipErrorNotMappedAsPointer = 213
        ENUMERATOR :: hipErrorECCNotCorrectable = 214
        ENUMERATOR :: hipErrorUnsupportedLimit = 215
        ENUMERATOR :: hipErrorContextAlreadyInUse = 216
        ENUMERATOR :: hipErrorPeerAccessUnsupported = 217
        ENUMERATOR :: hipErrorInvalidKernelFile = 218
        ENUMERATOR :: hipErrorInvalidGraphicsContext = 219
        ENUMERATOR :: hipErrorInvalidSource = 300
        ENUMERATOR :: hipErrorFileNotFound = 301
        ENUMERATOR :: hipErrorSharedObjectSymbolNotFound = 302
        ENUMERATOR :: hipErrorSharedObjectInitFailed = 303
        ENUMERATOR :: hipErrorOperatingSystem = 304
        ENUMERATOR :: hipErrorInvalidHandle = 400
        ENUMERATOR :: hipErrorInvalidResourceHandle = 400
        ENUMERATOR :: hipErrorIllegalState = 401
        ENUMERATOR :: hipErrorNotFound = 500
        ENUMERATOR :: hipErrorNotReady = 600
        ENUMERATOR :: hipErrorIllegalAddress = 700
        ENUMERATOR :: hipErrorLaunchOutOfResources = 701
        ENUMERATOR :: hipErrorLaunchTimeOut = 702
        ENUMERATOR :: hipErrorPeerAccessAlreadyEnabled = 704
        ENUMERATOR :: hipErrorPeerAccessNotEnabled = 705
        ENUMERATOR :: hipErrorSetOnActiveProcess = 708
        ENUMERATOR :: hipErrorContextIsDestroyed = 709
        ENUMERATOR :: hipErrorAssert = 710
        ENUMERATOR :: hipErrorHostMemoryAlreadyRegistered = 712
        ENUMERATOR :: hipErrorHostMemoryNotRegistered = 713
        ENUMERATOR :: hipErrorLaunchFailure = 719
        ENUMERATOR :: hipErrorCooperativeLaunchTooLarge = 720
        ENUMERATOR :: hipErrorNotSupported = 801
        ENUMERATOR :: hipErrorStreamCaptureUnsupported = 900
        ENUMERATOR :: hipErrorStreamCaptureInvalidated = 901
        ENUMERATOR :: hipErrorStreamCaptureMerge = 902
        ENUMERATOR :: hipErrorStreamCaptureUnmatched = 903
        ENUMERATOR :: hipErrorStreamCaptureUnjoined = 904
        ENUMERATOR :: hipErrorStreamCaptureIsolation = 905
        ENUMERATOR :: hipErrorStreamCaptureImplicit = 906
        ENUMERATOR :: hipErrorCapturedEvent = 907
        ENUMERATOR :: hipErrorStreamCaptureWrongThread = 908
        ENUMERATOR :: hipErrorGraphExecUpdateFailure = 910
        ENUMERATOR :: hipErrorUnknown = 999
        ENUMERATOR :: hipErrorRuntimeMemory = 1052
        ENUMERATOR :: hipErrorRuntimeOther = 1053
        ENUMERATOR :: hipErrorTbd
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipDeviceAttributeCudaCompatibleBegin = 0
        ENUMERATOR :: hipDeviceAttributeEccEnabled = hipDeviceAttributeCudaCompatibleBegin
        ENUMERATOR :: hipDeviceAttributeAccessPolicyMaxWindowSize
        ENUMERATOR :: hipDeviceAttributeAsyncEngineCount
        ENUMERATOR :: hipDeviceAttributeCanMapHostMemory
        ENUMERATOR :: hipDeviceAttributeCanUseHostPointerForRegisteredMem
        ENUMERATOR :: hipDeviceAttributeClockRate
        ENUMERATOR :: hipDeviceAttributeComputeMode
        ENUMERATOR :: hipDeviceAttributeComputePreemptionSupported
        ENUMERATOR :: hipDeviceAttributeConcurrentKernels
        ENUMERATOR :: hipDeviceAttributeConcurrentManagedAccess
        ENUMERATOR :: hipDeviceAttributeCooperativeLaunch
        ENUMERATOR :: hipDeviceAttributeCooperativeMultiDeviceLaunch
        ENUMERATOR :: hipDeviceAttributeDeviceOverlap
        ENUMERATOR :: hipDeviceAttributeDirectManagedMemAccessFromHost
        ENUMERATOR :: hipDeviceAttributeGlobalL1CacheSupported
        ENUMERATOR :: hipDeviceAttributeHostNativeAtomicSupported
        ENUMERATOR :: hipDeviceAttributeIntegrated
        ENUMERATOR :: hipDeviceAttributeIsMultiGpuBoard
        ENUMERATOR :: hipDeviceAttributeKernelExecTimeout
        ENUMERATOR :: hipDeviceAttributeL2CacheSize
        ENUMERATOR :: hipDeviceAttributeLocalL1CacheSupported
        ENUMERATOR :: hipDeviceAttributeLuid
        ENUMERATOR :: hipDeviceAttributeLuidDeviceNodeMask
        ENUMERATOR :: hipDeviceAttributeComputeCapabilityMajor
        ENUMERATOR :: hipDeviceAttributeManagedMemory
        ENUMERATOR :: hipDeviceAttributeMaxBlocksPerMultiProcessor
        ENUMERATOR :: hipDeviceAttributeMaxBlockDimX
        ENUMERATOR :: hipDeviceAttributeMaxBlockDimY
        ENUMERATOR :: hipDeviceAttributeMaxBlockDimZ
        ENUMERATOR :: hipDeviceAttributeMaxGridDimX
        ENUMERATOR :: hipDeviceAttributeMaxGridDimY
        ENUMERATOR :: hipDeviceAttributeMaxGridDimZ
        ENUMERATOR :: hipDeviceAttributeMaxSurface1D
        ENUMERATOR :: hipDeviceAttributeMaxSurface1DLayered
        ENUMERATOR :: hipDeviceAttributeMaxSurface2D
        ENUMERATOR :: hipDeviceAttributeMaxSurface2DLayered
        ENUMERATOR :: hipDeviceAttributeMaxSurface3D
        ENUMERATOR :: hipDeviceAttributeMaxSurfaceCubemap
        ENUMERATOR :: hipDeviceAttributeMaxSurfaceCubemapLayered
        ENUMERATOR :: hipDeviceAttributeMaxTexture1DWidth
        ENUMERATOR :: hipDeviceAttributeMaxTexture1DLayered
        ENUMERATOR :: hipDeviceAttributeMaxTexture1DLinear
        ENUMERATOR :: hipDeviceAttributeMaxTexture1DMipmap
        ENUMERATOR :: hipDeviceAttributeMaxTexture2DWidth
        ENUMERATOR :: hipDeviceAttributeMaxTexture2DHeight
        ENUMERATOR :: hipDeviceAttributeMaxTexture2DGather
        ENUMERATOR :: hipDeviceAttributeMaxTexture2DLayered
        ENUMERATOR :: hipDeviceAttributeMaxTexture2DLinear
        ENUMERATOR :: hipDeviceAttributeMaxTexture2DMipmap
        ENUMERATOR :: hipDeviceAttributeMaxTexture3DWidth
        ENUMERATOR :: hipDeviceAttributeMaxTexture3DHeight
        ENUMERATOR :: hipDeviceAttributeMaxTexture3DDepth
        ENUMERATOR :: hipDeviceAttributeMaxTexture3DAlt
        ENUMERATOR :: hipDeviceAttributeMaxTextureCubemap
        ENUMERATOR :: hipDeviceAttributeMaxTextureCubemapLayered
        ENUMERATOR :: hipDeviceAttributeMaxThreadsDim
        ENUMERATOR :: hipDeviceAttributeMaxThreadsPerBlock
        ENUMERATOR :: hipDeviceAttributeMaxThreadsPerMultiProcessor
        ENUMERATOR :: hipDeviceAttributeMaxPitch
        ENUMERATOR :: hipDeviceAttributeMemoryBusWidth
        ENUMERATOR :: hipDeviceAttributeMemoryClockRate
        ENUMERATOR :: hipDeviceAttributeComputeCapabilityMinor
        ENUMERATOR :: hipDeviceAttributeMultiGpuBoardGroupID
        ENUMERATOR :: hipDeviceAttributeMultiprocessorCount
        ENUMERATOR :: hipDeviceAttributeName
        ENUMERATOR :: hipDeviceAttributePageableMemoryAccess
        ENUMERATOR :: hipDeviceAttributePageableMemoryAccessUsesHostPageTables
        ENUMERATOR :: hipDeviceAttributePciBusId
        ENUMERATOR :: hipDeviceAttributePciDeviceId
        ENUMERATOR :: hipDeviceAttributePciDomainID
        ENUMERATOR :: hipDeviceAttributePersistingL2CacheMaxSize
        ENUMERATOR :: hipDeviceAttributeMaxRegistersPerBlock
        ENUMERATOR :: hipDeviceAttributeMaxRegistersPerMultiprocessor
        ENUMERATOR :: hipDeviceAttributeReservedSharedMemPerBlock
        ENUMERATOR :: hipDeviceAttributeMaxSharedMemoryPerBlock
        ENUMERATOR :: hipDeviceAttributeSharedMemPerBlockOptin
        ENUMERATOR :: hipDeviceAttributeSharedMemPerMultiprocessor
        ENUMERATOR :: hipDeviceAttributeSingleToDoublePrecisionPerfRatio
        ENUMERATOR :: hipDeviceAttributeStreamPrioritiesSupported
        ENUMERATOR :: hipDeviceAttributeSurfaceAlignment
        ENUMERATOR :: hipDeviceAttributeTccDriver
        ENUMERATOR :: hipDeviceAttributeTextureAlignment
        ENUMERATOR :: hipDeviceAttributeTexturePitchAlignment
        ENUMERATOR :: hipDeviceAttributeTotalConstantMemory
        ENUMERATOR :: hipDeviceAttributeTotalGlobalMem
        ENUMERATOR :: hipDeviceAttributeUnifiedAddressing
        ENUMERATOR :: hipDeviceAttributeUuid
        ENUMERATOR :: hipDeviceAttributeWarpSize
        ENUMERATOR :: hipDeviceAttributeCudaCompatibleEND = 9999
        ENUMERATOR :: hipDeviceAttributeAmdSpecificBegin = 10000
        ENUMERATOR :: hipDeviceAttributeClockInstructionRate = hipDeviceAttributeAmdSpecificBegin
        ENUMERATOR :: hipDeviceAttributeArch
        ENUMERATOR :: hipDeviceAttributeMaxSharedMemoryPerMultiprocessor
        ENUMERATOR :: hipDeviceAttributeGcnArch
        ENUMERATOR :: hipDeviceAttributeGcnArchName
        ENUMERATOR :: hipDeviceAttributeHdpMemFlushCntl
        ENUMERATOR :: hipDeviceAttributeHdpRegFlushCntl
        ENUMERATOR :: hipDeviceAttributeCooperativeMultiDeviceUnmatchedFunc
        ENUMERATOR :: hipDeviceAttributeCooperativeMultiDeviceUnmatchedGridDim
        ENUMERATOR :: hipDeviceAttributeCooperativeMultiDeviceUnmatchedBlockDim
        ENUMERATOR :: hipDeviceAttributeCooperativeMultiDeviceUnmatchedSharedMem
        ENUMERATOR :: hipDeviceAttributeIsLargeBar
        ENUMERATOR :: hipDeviceAttributeAsicRevision
        ENUMERATOR :: hipDeviceAttributeCanUseStreamWaitValue
        ENUMERATOR :: hipDeviceAttributeImageSupport
        ENUMERATOR :: hipDeviceAttributeAmdSpecificEND = 19999
        ENUMERATOR :: hipDeviceAttributeVENDorSpecificBegin = 20000
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipComputeModeDefault = 0
        ENUMERATOR :: hipComputeModeExclusive = 1
        ENUMERATOR :: hipComputeModeProhibited = 2
        ENUMERATOR :: hipComputeModeExclusiveProcess = 3
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipDevP2PAttrPerformanceRank = 0
        ENUMERATOR :: hipDevP2PAttrAccessSupported
        ENUMERATOR :: hipDevP2PAttrNativeAtomicSupported
        ENUMERATOR :: hipDevP2PAttrHipArrayAccessSupported
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipLimitPrintfFifoSize = 1
        ENUMERATOR :: hipLimitMallocHeapSize = 2
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipMemAdviseSetReadMostly = 1
        ENUMERATOR :: hipMemAdviseUnsetReadMostly = 2
        ENUMERATOR :: hipMemAdviseSetPreferredLocation = 3
        ENUMERATOR :: hipMemAdviseUnsetPreferredLocation = 4
        ENUMERATOR :: hipMemAdviseSetAccessedBy = 5
        ENUMERATOR :: hipMemAdviseUnsetAccessedBy = 6
        ENUMERATOR :: hipMemAdviseSetCoarseGrain = 100
        ENUMERATOR :: hipMemAdviseUnsetCoarseGrain = 101
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipMemRangeCoherencyModeFineGrain = 0
        ENUMERATOR :: hipMemRangeCoherencyModeCoarseGrain = 1
        ENUMERATOR :: hipMemRangeCoherencyModeIndeterminate = 2
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipMemRangeAttributeReadMostly = 1
        ENUMERATOR :: hipMemRangeAttributePreferredLocation = 2
        ENUMERATOR :: hipMemRangeAttributeAccessedBy = 3
        ENUMERATOR :: hipMemRangeAttributeLastPrefetchLocation = 4
        ENUMERATOR :: hipMemRangeAttributeCoherencyMode = 100
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipJitOptionMaxRegisters = 0
        ENUMERATOR :: hipJitOptionThreadsPerBlock
        ENUMERATOR :: hipJitOptionWallTime
        ENUMERATOR :: hipJitOptionInfoLogBuffer
        ENUMERATOR :: hipJitOptionInfoLogBufferSizeBytes
        ENUMERATOR :: hipJitOptionErrorLogBuffer
        ENUMERATOR :: hipJitOptionErrorLogBufferSizeBytes
        ENUMERATOR :: hipJitOptionOptimizationLevel
        ENUMERATOR :: hipJitOptionTargetFromContext
        ENUMERATOR :: hipJitOptionTarget
        ENUMERATOR :: hipJitOptionFallbackStrategy
        ENUMERATOR :: hipJitOptionGenerateDebugInfo
        ENUMERATOR :: hipJitOptionLogVerbose
        ENUMERATOR :: hipJitOptionGenerateLineInfo
        ENUMERATOR :: hipJitOptionCacheMode
        ENUMERATOR :: hipJitOptionSm3xOpt
        ENUMERATOR :: hipJitOptionFastCompile
        ENUMERATOR :: hipJitOptionNumOptions
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipFuncAttributeMaxDynamicSharedMemorySize = 8
        ENUMERATOR :: hipFuncAttributePreferredSharedMemoryCarveout = 9
        ENUMERATOR :: hipFuncAttributeMax
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipFuncCachePreferNone
        ENUMERATOR :: hipFuncCachePreferShared
        ENUMERATOR :: hipFuncCachePreferL1
        ENUMERATOR :: hipFuncCachePreferEqual
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipSharedMemBankSizeDefault
        ENUMERATOR :: hipSharedMemBankSizeFourByte
        ENUMERATOR :: hipSharedMemBankSizeEightByte
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipExternalMemoryHandleTYPEOpaqueFd = 1
        ENUMERATOR :: hipExternalMemoryHandleTYPEOpaqueWin32 = 2
        ENUMERATOR :: hipExternalMemoryHandleTYPEOpaqueWin32Kmt = 3
        ENUMERATOR :: hipExternalMemoryHandleTYPED3D12Heap = 4
        ENUMERATOR :: hipExternalMemoryHandleTYPED3D12Resource = 5
        ENUMERATOR :: hipExternalMemoryHandleTYPED3D11Resource = 6
        ENUMERATOR :: hipExternalMemoryHandleTYPED3D11ResourceKmt = 7
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipExternalSemaphoreHandleTYPEOpaqueFd = 1
        ENUMERATOR :: hipExternalSemaphoreHandleTYPEOpaqueWin32 = 2
        ENUMERATOR :: hipExternalSemaphoreHandleTYPEOpaqueWin32Kmt = 3
        ENUMERATOR :: hipExternalSemaphoreHandleTYPED3D12Fence = 4
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipGLDeviceListAll = 1
        ENUMERATOR :: hipGLDeviceListCurrentFrame = 2
        ENUMERATOR :: hipGLDeviceListNextFrame = 3
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipGraphicsRegisterFlagsNone = 0
        ENUMERATOR :: hipGraphicsRegisterFlagsReadOnly = 1
        ENUMERATOR :: hipGraphicsRegisterFlagsWriteDiscard = 2
        ENUMERATOR :: hipGraphicsRegisterFlagsSurfaceLoadStore = 4
        ENUMERATOR :: hipGraphicsRegisterFlagsTextureGather = 8
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipGraphNodeTYPEKernel = 1
        ENUMERATOR :: hipGraphNodeTYPEMemcpy = 2
        ENUMERATOR :: hipGraphNodeTYPEMemset = 3
        ENUMERATOR :: hipGraphNodeTYPEHost = 4
        ENUMERATOR :: hipGraphNodeTYPEGraph = 5
        ENUMERATOR :: hipGraphNodeTYPEEmpty = 6
        ENUMERATOR :: hipGraphNodeTYPEWaitEvent = 7
        ENUMERATOR :: hipGraphNodeTYPEEventRecord = 8
        ENUMERATOR :: hipGraphNodeTYPEMemcpy1D = 9
        ENUMERATOR :: hipGraphNodeTYPEMemcpyFromSymbol = 10
        ENUMERATOR :: hipGraphNodeTYPEMemcpyToSymbol = 11
        ENUMERATOR :: hipGraphNodeTYPECount
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipGraphExecUpdateSuccess = 0
        ENUMERATOR :: hipGraphExecUpdateError = 1
        ENUMERATOR :: hipGraphExecUpdateErrorTopologyChanged = 2
        ENUMERATOR :: hipGraphExecUpdateErrorNodeTYPEChanged = 3
        ENUMERATOR :: hipGraphExecUpdateErrorFUNCTIONChanged = 4
        ENUMERATOR :: hipGraphExecUpdateErrorParametersChanged = 5
        ENUMERATOR :: hipGraphExecUpdateErrorNotSupported = 6
        ENUMERATOR :: hipGraphExecUpdateErrorUnsupportedFUNCTIONChange = 7
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipStreamCaptureModeGlobal = 0
        ENUMERATOR :: hipStreamCaptureModeThreadLocal
        ENUMERATOR :: hipStreamCaptureModeRelaxed
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipStreamCaptureStatusNone = 0
        ENUMERATOR :: hipStreamCaptureStatusActive
        ENUMERATOR :: hipStreamCaptureStatusInvalidated
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipStreamAddCaptureDepENDencies = 0
        ENUMERATOR :: hipStreamSetCaptureDepENDencies
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipChannelFormatKindSigned = 0
        ENUMERATOR :: hipChannelFormatKindUnsigned = 1
        ENUMERATOR :: hipChannelFormatKindFloat = 2
        ENUMERATOR :: hipChannelFormatKindNone = 3
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: HIP_AD_FORMAT_UNSIGNED_INT8 = 1
        ENUMERATOR :: HIP_AD_FORMAT_UNSIGNED_INT16 = 2
        ENUMERATOR :: HIP_AD_FORMAT_UNSIGNED_INT32 = 3
        ENUMERATOR :: HIP_AD_FORMAT_SIGNED_INT8 = 8
        ENUMERATOR :: HIP_AD_FORMAT_SIGNED_INT16 = 9
        ENUMERATOR :: HIP_AD_FORMAT_SIGNED_INT32 = 10
        ENUMERATOR :: HIP_AD_FORMAT_HALF = 16
        ENUMERATOR :: HIP_AD_FORMAT_FLOAT = 32
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipResourceTYPEArray = 0
        ENUMERATOR :: hipResourceTYPEMipmappedArray = 1
        ENUMERATOR :: hipResourceTYPELinear = 2
        ENUMERATOR :: hipResourceTYPEPitch2D = 3
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: HIP_RESOURCE_TYPE_ARRAY = 0
        ENUMERATOR :: HIP_RESOURCE_TYPE_MIPMAPPED_ARRAY = 1
        ENUMERATOR :: HIP_RESOURCE_TYPE_LINEAR = 2
        ENUMERATOR :: HIP_RESOURCE_TYPE_PITCH2D = 3
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: HIP_TR_ADDRESS_MODE_WRAP = 0
        ENUMERATOR :: HIP_TR_ADDRESS_MODE_CLAMP = 1
        ENUMERATOR :: HIP_TR_ADDRESS_MODE_MIRROR = 2
        ENUMERATOR :: HIP_TR_ADDRESS_MODE_BORDER = 3
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: HIP_TR_FILTER_MODE_POINT = 0
        ENUMERATOR :: HIP_TR_FILTER_MODE_LINEAR = 1
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipResViewFormatNone = 0
        ENUMERATOR :: hipResViewFormatUnsignedChar1 = 1
        ENUMERATOR :: hipResViewFormatUnsignedChar2 = 2
        ENUMERATOR :: hipResViewFormatUnsignedChar4 = 3
        ENUMERATOR :: hipResViewFormatSignedChar1 = 4
        ENUMERATOR :: hipResViewFormatSignedChar2 = 5
        ENUMERATOR :: hipResViewFormatSignedChar4 = 6
        ENUMERATOR :: hipResViewFormatUnsignedShort1 = 7
        ENUMERATOR :: hipResViewFormatUnsignedShort2 = 8
        ENUMERATOR :: hipResViewFormatUnsignedShort4 = 9
        ENUMERATOR :: hipResViewFormatSignedShort1 = 10
        ENUMERATOR :: hipResViewFormatSignedShort2 = 11
        ENUMERATOR :: hipResViewFormatSignedShort4 = 12
        ENUMERATOR :: hipResViewFormatUnsignedInt1 = 13
        ENUMERATOR :: hipResViewFormatUnsignedInt2 = 14
        ENUMERATOR :: hipResViewFormatUnsignedInt4 = 15
        ENUMERATOR :: hipResViewFormatSignedInt1 = 16
        ENUMERATOR :: hipResViewFormatSignedInt2 = 17
        ENUMERATOR :: hipResViewFormatSignedInt4 = 18
        ENUMERATOR :: hipResViewFormatHalf1 = 19
        ENUMERATOR :: hipResViewFormatHalf2 = 20
        ENUMERATOR :: hipResViewFormatHalf4 = 21
        ENUMERATOR :: hipResViewFormatFloat1 = 22
        ENUMERATOR :: hipResViewFormatFloat2 = 23
        ENUMERATOR :: hipResViewFormatFloat4 = 24
        ENUMERATOR :: hipResViewFormatUnsignedBlockCompressed1 = 25
        ENUMERATOR :: hipResViewFormatUnsignedBlockCompressed2 = 26
        ENUMERATOR :: hipResViewFormatUnsignedBlockCompressed3 = 27
        ENUMERATOR :: hipResViewFormatUnsignedBlockCompressed4 = 28
        ENUMERATOR :: hipResViewFormatSignedBlockCompressed4 = 29
        ENUMERATOR :: hipResViewFormatUnsignedBlockCompressed5 = 30
        ENUMERATOR :: hipResViewFormatSignedBlockCompressed5 = 31
        ENUMERATOR :: hipResViewFormatUnsignedBlockCompressed6H = 32
        ENUMERATOR :: hipResViewFormatSignedBlockCompressed6H = 33
        ENUMERATOR :: hipResViewFormatUnsignedBlockCompressed7 = 34
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_NONE = 0
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UINT_1X8 = 1
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UINT_2X8 = 2
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UINT_4X8 = 3
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SINT_1X8 = 4
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SINT_2X8 = 5
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SINT_4X8 = 6
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UINT_1X16 = 7
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UINT_2X16 = 8
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UINT_4X16 = 9
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SINT_1X16 = 10
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SINT_2X16 = 11
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SINT_4X16 = 12
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UINT_1X32 = 13
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UINT_2X32 = 14
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UINT_4X32 = 15
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SINT_1X32 = 16
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SINT_2X32 = 17
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SINT_4X32 = 18
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_FLOAT_1X16 = 19
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_FLOAT_2X16 = 20
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_FLOAT_4X16 = 21
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_FLOAT_1X32 = 22
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_FLOAT_2X32 = 23
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_FLOAT_4X32 = 24
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UNSIGNED_BC1 = 25
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UNSIGNED_BC2 = 26
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UNSIGNED_BC3 = 27
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UNSIGNED_BC4 = 28
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SIGNED_BC4 = 29
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UNSIGNED_BC5 = 30
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SIGNED_BC5 = 31
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UNSIGNED_BC6H = 32
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_SIGNED_BC6H = 33
        ENUMERATOR :: HIP_RES_VIEW_FORMAT_UNSIGNED_BC7 = 34
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipMemcpyHostToHost = 0
        ENUMERATOR :: hipMemcpyHostToDevice = 1
        ENUMERATOR :: hipMemcpyDeviceToHost = 2
        ENUMERATOR :: hipMemcpyDeviceToDevice = 3
        ENUMERATOR :: hipMemcpyDefault = 4
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: HIP_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK
        ENUMERATOR :: HIP_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES
        ENUMERATOR :: HIP_FUNC_ATTRIBUTE_CONST_SIZE_BYTES
        ENUMERATOR :: HIP_FUNC_ATTRIBUTE_LOCAL_SIZE_BYTES
        ENUMERATOR :: HIP_FUNC_ATTRIBUTE_NUM_REGS
        ENUMERATOR :: HIP_FUNC_ATTRIBUTE_PTX_VERSION
        ENUMERATOR :: HIP_FUNC_ATTRIBUTE_BINARY_VERSION
        ENUMERATOR :: HIP_FUNC_ATTRIBUTE_CACHE_MODE_CA
        ENUMERATOR :: HIP_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES
        ENUMERATOR :: HIP_FUNC_ATTRIBUTE_PREFERRED_SHARED_MEMORY_CARVEOUT
        ENUMERATOR :: HIP_FUNC_ATTRIBUTE_MAX
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_CONTEXT = 1
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_MEMORY_TYPE
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_DEVICE_POINTER
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_HOST_POINTER
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_P2P_TOKENS
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_SYNC_MEMOPS
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_BUFFER_ID
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_IS_MANAGED
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_DEVICE_ORDINAL
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_IS_LEGACY_HIP_IPC_CAPABLE
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_RANGE_START_ADDR
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_RANGE_SIZE
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_MAPPED
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_ALLOWED_HANDLE_TYPES
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_IS_GPU_DIRECT_RDMA_CAPABLE
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_ACCESS_FLAGS
        ENUMERATOR :: HIP_POINTER_ATTRIBUTE_MEMPOOL_HANDLE
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipAddressModeWrap = 0
        ENUMERATOR :: hipAddressModeClamp = 1
        ENUMERATOR :: hipAddressModeMirror = 2
        ENUMERATOR :: hipAddressModeBorder = 3
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipFilterModePoint = 0
        ENUMERATOR :: hipFilterModeLinear = 1
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: hipReadModeElementTYPE = 0
        ENUMERATOR :: hipReadModeNormalizedFloat = 1
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: HIP_R_16F = 2
        ENUMERATOR :: HIP_R_32F = 0
        ENUMERATOR :: HIP_R_64F = 1
        ENUMERATOR :: HIP_C_16F = 6
        ENUMERATOR :: HIP_C_32F = 4
        ENUMERATOR :: HIP_C_64F = 5
      END ENUM
    
      ENUM,BIND(c)
        ENUMERATOR :: HIP_LIBRARY_MAJOR_VERSION
        ENUMERATOR :: HIP_LIBRARY_MINOR_VERSION
        ENUMERATOR :: HIP_LIBRARY_PATCH_LEVEL
      END ENUM
      
    END MODULE SELF_HIP_enums
    