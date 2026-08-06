! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2026 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

program test
  !! Guard test (registered with WILL_FAIL): Read_HOHQMesh must abort
  !! with a nonzero stop code when the mesh file cannot be opened.

  use SELF_Mesh_3D

  implicit none

  type(Mesh3D),target :: mesh

  call mesh%Read_HOHQMesh("this_mesh_file_does_not_exist.mesh")

  ! Unreachable when the guard works; exiting cleanly here would make
  ! the WILL_FAIL test fail.
  print*,"FAIL: Read_HOHQMesh did not abort on a missing file"

endprogram test
