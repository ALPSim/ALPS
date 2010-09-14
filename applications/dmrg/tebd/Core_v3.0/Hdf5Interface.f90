!    Copyright (C) 2010  M. L. Wall and L. D. Carr, Colorado School of Mines
!    This file is part of OpenSourceTEBD.
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
MODULE Hdf5Interface
!
! Purpose: Interface to the hdf5 fortran API
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!       8/18/10  M. L. Wall	alpha release
!       9/13/10  M. L. Wall	Changed quantum number i/o
!
USE GlobalData
USE LinearOps
USE HamiOps
USE StateOps
USE LocalOps
USE ObOps
USE HDF5
IMPLICIT NONE
INTEGER :: error
INTEGER(HID_T) :: outputfileId, timeStepsGrpID, paramsgrpId

INTERFACE WriteDataToGroup
	MODULE PROCEDURE WriteDataToGroup_I, WriteDataToGroup_Iv, WriteDataToGroup_Im, WriteDataToGroup_r,&
	 WriteDataToGroup_rv, WriteDataToGroup_rm, WriteDataToGroup_c,  WriteDataToGroup_L, WriteDataToGroup_cv
END INTERFACE WriteDataToGroup

INTERFACE WriteDataToMeanGroup
	MODULE PROCEDURE WriteDataToMeanGroup_I, WriteDataToMeanGroup_Iv, WriteDataToMeanGroup_Im, WriteDataToMeanGroup_r,&
	 WriteDataToMeanGroup_rv, WriteDataToMeanGroup_rm, WriteDataToMeanGroup_c,  WriteDataToMeanGroup_L, WriteDataToMeanGroup_cv
END INTERFACE WriteDataToMeanGroup

CONTAINS 

SUBROUTINE InitializeHdf5()
IMPLICIT NONE

 CALL h5open_f(error)!Initialize FORTRAN interface.

END SUBROUTINE InitializeHdf5

SUBROUTINE FinalizeHdf5()
IMPLICIT NONE

 CALL h5close_f(error) !    Close FORTRAN interface.

END SUBROUTINE FinalizeHdf5

SUBROUTINE OpenHdf5File(fileName,fileId, openKind)
!
!Purpose : Open an hdf5 File.  
!OpenKind='rw' speficies readwrite, otherwise the file is deleted if it already exists.
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) :: fileName
CHARACTER(len=*), INTENT(IN), OPTIONAL :: openKind
INTEGER(HID_T), INTENT(INOUT) :: fileId

IF(PRESENT(openKind)) THEN
	IF(openKind=='rw') THEN
		CALL h5fopen_f(filename, H5F_ACC_RDWR_F, fileid, error)
	ELSE
		PRINT *, 'Unknown openKind option in OpenHdf5File!'
	END IF
ELSE
	CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, fileid, error)
END IF

END SUBROUTINE OpenHdf5File

SUBROUTINE CloseHdf5File(fileId)
!
!Purpose : Terminate access to and close an hdf5 File.  
!
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: fileId

 CALL h5fclose_f(fileid, error) ! Terminate access to the file.

END SUBROUTINE CloseHdf5File

SUBROUTINE CreateGroup(groupName, groupId, fileId)
!
!Purpose: Create a group called groupName in the file labelled by fileId
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: fileId
CHARACTER(len=*), INTENT(IN) :: groupName
INTEGER, INTENT(INOUT) :: groupId

 CALL h5gcreate_f(fileid, groupname, groupid, error) ! Create a group named "groupname" in the file.

END SUBROUTINE CreateGroup

SUBROUTINE OpenGroup(groupName, groupId, fileId)
!
!Purpose: Create a group called groupName in the file labelled by fileId
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: fileId
CHARACTER(len=*), INTENT(IN) :: groupName
INTEGER, INTENT(INOUT) :: groupId

 CALL h5gopen_f(fileId,groupName,groupid, error) ! Create a group named "groupname" in the file.

END SUBROUTINE OpenGroup

SUBROUTINE CreateSubGroup(subGroupName, subGroupId, groupId)
!
!Purpose: create a subgroup inside of an existing group
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) :: subGroupName
INTEGER, INTENT(IN) :: groupId
INTEGER, INTENT(INOUT) :: subGroupId

 CALL h5gcreate_f(groupid, subgroupName, subgroupid, error) ! Create group "subgroupName1" in group "groupname" using relative name.

END SUBROUTINE CreateSubGroup

SUBROUTINE WriteDataToGroup_I(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) :: dataspaceID, dataSetID
INTEGER, INTENT(IN) :: myData
INTEGER(HSIZE_T) :: dims(1)

 dims=1
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL h5dcreate_f(groupId, dataName, H5T_NATIVE_INTEGER, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_INTEGER, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.

END SUBROUTINE WriteDataToGroup_I

SUBROUTINE WriteDataToGroup_Iv(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) ::  dataspaceID, dataSetID
INTEGER, INTENT(IN) :: myData(:)
INTEGER(HSIZE_T) :: dims(1)

 dims=SIZE(myData)
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL h5dcreate_f(groupId, dataName, H5T_NATIVE_INTEGER, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_INTEGER, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.

END SUBROUTINE WriteDataToGroup_Iv

SUBROUTINE WriteDataToGroup_Im(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) ::  dataspaceID, dataSetID
INTEGER, INTENT(IN) :: myData(:,:)
INTEGER(HSIZE_T) :: dims(2)

 dims=(/SIZE(myData,1), SIZE(myData,2)/)
 rank=2
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL h5dcreate_f(groupId, dataName, H5T_NATIVE_INTEGER, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_INTEGER, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.

END SUBROUTINE WriteDataToGroup_Im

SUBROUTINE WriteDataToGroup_r(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) :: dataspaceID, dataSetID
REAL(KIND=8), INTENT(IN) :: myData
INTEGER(HSIZE_T) :: dims(1)

 dims=1
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL h5dcreate_f(groupId, dataName, H5T_NATIVE_DOUBLE, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_DOUBLE, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.

END SUBROUTINE WriteDataToGroup_r

SUBROUTINE WriteDataToGroup_rv(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) ::  dataspaceID, dataSetID
REAL(KIND=8), INTENT(IN) :: myData(:)
INTEGER(HSIZE_T) :: dims(1)

 dims=SIZE(myData)
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL h5dcreate_f(groupId, dataName, H5T_NATIVE_DOUBLE, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_DOUBLE, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.

END SUBROUTINE WriteDataToGroup_rv

SUBROUTINE WriteDataToGroup_rm(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) ::  dataspaceID, dataSetID
REAL(KIND=8), INTENT(IN) :: myData(:,:)
INTEGER(HSIZE_T) :: dims(2)

 dims=(/SIZE(myData,1), SIZE(myData,2)/)
 rank=2
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL h5dcreate_f(groupId, dataName, H5T_NATIVE_DOUBLE, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_DOUBLE, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.

END SUBROUTINE WriteDataToGroup_rm

SUBROUTINE WriteDataToGroup_c(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) :: dataspaceID, dataSetID
CHARACTER(len=*), INTENT(IN) :: myData
INTEGER(HSIZE_T) :: dims(1)
INTEGER(HID_T) :: strtype

 CALL h5tcopy_f(H5T_NATIVE_CHARACTER,strtype,error)
 CALL h5tset_size_f(strtype,LEN_TRIM(myData),error)
 dims=1
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL h5dcreate_f(groupId, dataName, strtype, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, strtype, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.

END SUBROUTINE WriteDataToGroup_c

SUBROUTINE WriteDataToGroup_cv(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) :: dataspaceID, dataSetID
CHARACTER(len=*), INTENT(IN) :: myData(:)
INTEGER(HSIZE_T) :: dims(1)
INTEGER(HID_T) :: strtype

 CALL h5tcopy_f(H5T_NATIVE_CHARACTER,strtype,error)
 CALL h5tset_size_f(strtype,10,error)
 dims=SIZE(myData)
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL h5dcreate_f(groupId, dataName, strtype, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, strtype, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.

END SUBROUTINE WriteDataToGroup_cv

SUBROUTINE WriteDataToGroup_L(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
INTEGER(HID_T), INTENT(IN) :: groupId
LOGICAL, INTENT(IN) :: myData
CHARACTER(len=5) :: strData

IF(myData) THEN
	strDAta='TRUE'
ELSE
	strDAta='FALSE'
END IF

CALL WriteDataToGroup(strData, groupId, dataName)

END SUBROUTINE WriteDataToGroup_L

SUBROUTINE WriteDataToMeanGroup_I(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) :: dataspaceID, dataSetID, meangrpID
INTEGER, INTENT(IN) :: myData
INTEGER(HSIZE_T) :: dims(1)

 dims=1
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL CreateSubGroup("mean",meangrpID,groupId)
 CALL h5dcreate_f(meangrpID, dataName, H5T_NATIVE_INTEGER, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_INTEGER, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.
 CALL CloseGroup(meangrpID)

END SUBROUTINE WriteDataToMeanGroup_I

SUBROUTINE WriteDataToMeanGroup_Iv(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) ::  dataspaceID, dataSetID, meangrpID
INTEGER, INTENT(IN) :: myData(:)
INTEGER(HSIZE_T) :: dims(1)

 dims=SIZE(myData)
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL CreateSubGroup("mean",meangrpID,groupId)
 CALL h5dcreate_f(meangrpID, dataName, H5T_NATIVE_INTEGER, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_INTEGER, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.
 CALL CloseGroup(meangrpID)

END SUBROUTINE WriteDataToMeanGroup_Iv

SUBROUTINE WriteDataToMeanGroup_Im(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) ::  dataspaceID, dataSetID, meangrpID
INTEGER, INTENT(IN) :: myData(:,:)
INTEGER(HSIZE_T) :: dims(2)

 dims=(/SIZE(myData,1), SIZE(myData,2)/)
 rank=2
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL CreateSubGroup("mean",meangrpID,groupId)
 CALL h5dcreate_f(meangrpID, dataName, H5T_NATIVE_INTEGER, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_INTEGER, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.
 CALL CloseGroup(meangrpID)

END SUBROUTINE WriteDataToMeanGroup_Im

SUBROUTINE WriteDataToMeanGroup_r(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) :: dataspaceID, dataSetID, meangrpID
REAL(KIND=8), INTENT(IN) :: myData
INTEGER(HSIZE_T) :: dims(1)

 dims=1
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL CreateSubGroup("mean",meangrpID,groupId)
 CALL h5dcreate_f(meangrpID, dataName, H5T_NATIVE_DOUBLE, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_DOUBLE, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.
 CALL CloseGroup(meangrpID)

END SUBROUTINE WriteDataToMeanGroup_r

SUBROUTINE WriteDataToMeanGroup_rv(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) ::  dataspaceID, dataSetID, meangrpID
REAL(KIND=8), INTENT(IN) :: myData(:)
INTEGER(HSIZE_T) :: dims(1)

 dims=SIZE(myData)
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL CreateSubGroup("mean",meangrpID,groupId)
 CALL h5dcreate_f(meangrpID, dataName, H5T_NATIVE_DOUBLE, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_DOUBLE, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.
 CALL CloseGroup(meangrpID)

END SUBROUTINE WriteDataToMeanGroup_rv

SUBROUTINE WriteDataToMeanGroup_rm(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) ::  dataspaceID, dataSetID, meangrpID
REAL(KIND=8), INTENT(IN) :: myData(:,:)
INTEGER(HSIZE_T) :: dims(2)

 dims=(/SIZE(myData,1), SIZE(myData,2)/)
 rank=2
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL CreateSubGroup("mean",meangrpID,groupId)
 CALL h5dcreate_f(meangrpID, dataName, H5T_NATIVE_DOUBLE, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, H5T_NATIVE_DOUBLE, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.
 CALL CloseGroup(meangrpID)

END SUBROUTINE WriteDataToMeanGroup_rm

SUBROUTINE WriteDataToMeanGroup_c(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) :: dataspaceID, dataSetID, meangrpID
CHARACTER(len=*), INTENT(IN) :: myData
INTEGER(HSIZE_T) :: dims(1)
INTEGER(HID_T) :: strtype

 CALL h5tcopy_f(H5T_NATIVE_CHARACTER,strtype,error)
 CALL h5tset_size_f(strtype,LEN_TRIM(myData),error)
 dims=1
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL CreateSubGroup("mean",meangrpID,groupId)
 CALL h5dcreate_f(meangrpID, dataName, strtype, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, strtype, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.
 CALL CloseGroup(meangrpID)

END SUBROUTINE WriteDataToMeanGroup_c


SUBROUTINE WriteDataToMeanGroup_cv(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
CHARACTER(len=132) :: dsetName
INTEGER(HID_T), INTENT(IN) :: groupId
INTEGER :: rank
INTEGER(HID_T) :: dataspaceID, dataSetID, meangrpID
CHARACTER(len=*), INTENT(IN) :: myData(:)
INTEGER(HSIZE_T) :: dims(1)
INTEGER(HID_T) :: strtype

 CALL h5tcopy_f(H5T_NATIVE_CHARACTER,strtype,error)
 CALL h5tset_size_f(strtype,10,error)
 dims=SIZE(myData)
 rank=1
 CALL h5screate_simple_f(rank, dims, dataspaceID, error) ! Create the data space for the first dataset
 CALL CreateSubGroup("mean",meangrpID,groupId)
 CALL h5dcreate_f(meangrpID, dataName, strtype, dataspaceID, dataSetID, error) ! Create a dataset in group "MyGroup" with default properties.
 CALL h5dwrite_f(dataSetID, strtype, myData, dims, error)
 CALL h5sclose_f(dataspaceID, error) ! Close the dataspace for the first dataset.
 CALL h5dclose_f(dataSetID, error) ! Close the first dataset.
 CALL CloseGroup(meangrpID)

END SUBROUTINE WriteDataToMeanGroup_cv

SUBROUTINE WriteDataToMeanGroup_L(myData, groupId, dataName)
!
!Purpose : Write myData to the group with absolute group name absGroupName and name dataStub in the file indicated by fileName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) ::  dataName
INTEGER(HID_T), INTENT(IN) :: groupId
LOGICAL, INTENT(IN) :: myData
CHARACTER(len=5) :: strData

IF(myData) THEN
	strDAta='TRUE'
ELSE
	strDAta='FALSE'
END IF

CALL WriteDataToMeanGroup(strData, groupId, dataName)

END SUBROUTINE WriteDataToMeanGroup_L

SUBROUTINE CloseGroup(groupId)
!
!Purpose:Close a group.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: groupId

 CALL h5gclose_f(groupid, error)! Close the group.

END SUBROUTINE CloseGroup


SUBROUTINE WriteSystemSettings()
!
!Purpose: Write the simulation parameters into a group in the hdf5 file
!
IMPLICIT NONE
REAL(KIND=rKind) :: locQ


CALL InitializeHdf5()
outputName=nmlName(1:LEN_TRIM(nmlName)-4)//".h5"
CALL OpenHdf5File(outputName,outputFileID)
CALL CreateGroup("parameters", paramsgrpId, outputFileID) !Create parameters group
CALL WriteDataToGroup(systemSize,paramsgrpId,"L")
CALL WriteDataToGroup(HamiType,paramsgrpId,"MODEL")
CALL WriteDataToGroup(initialState,paramsgrpId,"INITIAL STATE")
CALL WriteDataToGroup(rtp,paramsgrpId,"RTP")
CALL WriteDataToGroup(qSwitch,paramsgrpId,"QNUMBERS USED")
IF(qSwitch) THEN
	IF(Hamitype=='spin') THEN
		locQ=totQ-systemSize*spin
		CALL WriteDataToGroup(locQ,paramsgrpId,qType)
	ELSE
		CALL WriteDataToGroup(totQ,paramsgrpId,qType)
	END IF
END IF
CALL WriteDataToGroup(numThr,paramsgrpId,"NUM THREADS")
CALL WriteDataToGroup(trotterOrder,paramsgrpId,"TROTTER ORDER")
CALL WriteDataToGroup(chiLimit,paramsgrpId,"LIMITING RTP CHI VALUE")
CALL WriteDataToGroup(truncLimit,paramsgrpId,"LIMITING RTP TRUNC VALUE")
CALL WriteDataToGroup(simId,paramsgrpId,"SIMID")
CALL CloseGroup(paramsgrpID)
CALL CloseHdf5File(outputFileID)

END SUBROUTINE WriteSystemSettings


SUBROUTINE WriteITPSettings(Hamitype, Hparams)
!
!Purpose: Write the itp parameters into a group in the hdf5 file
!
IMPLICIT NONE
CHARACTER(len=*) :: HamiType
TYPE(HamiParams) :: Hparams

CALL OpenHdf5File(outputName,outputFileID, openKind='rw')
CALL OpenGroup("parameters", paramsgrpId, outputFileID) !Create parameters group
CALL WriteDataToGroup(numITP,paramsgrpId,"NUMITP")
CALL WriteDataToGroup(chivals,paramsgrpId,"CHIVALS")
CALL WriteDataToGroup(dtITPvals,paramsgrpId,"DTITPVALS")
CALL WriteDataToGroup(convCriterion,paramsgrpId,"CONVCRITERIA")

SELECT CASE(HamiType)
	CASE ('spin')
		CALL WriteDataToGroup(spin,paramsgrpId,"SPIN")
		CALL WriteDataToGroup(Hparams%sp%Jz,paramsgrpId,"ITP JZ")
		CALL WriteDataToGroup(Hparams%sp%Jxy,paramsgrpId,"ITP JXY")
		CALL WriteDataToGroup(Hparams%sp%h,paramsgrpId,"ITP H")
		CALL WriteDataToGroup(Hparams%sp%gam,paramsgrpId,"ITP GAMMA")
		CALL WriteDataToGroup(Hparams%sp%d,paramsgrpId,"ITP D")
		CALL WriteDataToGroup(Hparams%sp%k,paramsgrpId,"ITP K")
	CASE ('boson Hubbard')
		CALL WriteDataToGroup(Nmax,paramsgrpId,"Nmax")
		CALL WriteDataToGroup(Hparams%bp%t,paramsgrpId,"ITP T")
		CALL WriteDataToGroup(Hparams%bp%U,paramsgrpId,"ITP U")
		CALL WriteDataToGroup(Hparams%bp%V,paramsgrpId,"ITP V")
		CALL WriteDataToGroup(Hparams%bp%mu,paramsgrpId,"ITP MU")
	CASE ('hardcore boson')
		CALL WriteDataToGroup(Hparams%hcbp%t,paramsgrpId,"ITP T")
		CALL WriteDataToGroup(Hparams%hcbp%V,paramsgrpId,"ITP V")
		CALL WriteDataToGroup(Hparams%hcbp%mu,paramsgrpId,"ITP MU")
	CASE ('fermion Hubbard')
		CALL WriteDataToGroup(Hparams%fhp%t,paramsgrpId,"ITP T")
		CALL WriteDataToGroup(Hparams%fhp%U,paramsgrpId,"ITP U")
		CALL WriteDataToGroup(Hparams%fhp%V,paramsgrpId,"ITP V")
		CALL WriteDataToGroup(Hparams%fhp%mu,paramsgrpId,"ITP MU")
	CASE ('spinless fermions')
		CALL WriteDataToGroup(Hparams%sfp%t,paramsgrpId,"ITP T")
		CALL WriteDataToGroup(Hparams%sfp%V,paramsgrpId,"ITP V")
		CALL WriteDataToGroup(Hparams%sfp%mu,paramsgrpId,"ITP MU")
	CASE DEFAULT
		PRINT *, "Hamiltonian type not recognized in WriteITPSettings!"
		PRINT *, "Use 'spin', 'boson Hubbard', 'hardcore boson', 'fermion Hubbard', or 'spinless fermions'."
		STOP 
END SELECT
CALL CloseGroup(paramsgrpID)
CALL CloseHdf5File(outputFileID)


END SUBROUTINE WriteITPSettings

SUBROUTINE WriteRTPSettings(Hamitype,Hparams, rtpD)
!
!Purpose: Write the itp parameters into a group in the hdf5 file
!
IMPLICIT NONE
CHARACTER(len=*) :: HamiType
TYPE(rtpData) :: rtpD
TYPE(HamiParams) :: Hparams

CALL OpenHdf5File(outputName,outputFileID, openKind='rw')
CALL OpenGroup("parameters", paramsgrpId, outputFileID) !Open parameters group
SELECT CASE(HamiType)
	CASE ('spin')
		IF(.not.itp) THEN
			CALL WriteDataToGroup(spin,paramsgrpId,"SPIN")
		END IF
		CALL WriteDataToGroup(Hparams%sp%Jz,paramsgrpId,"JZ")
		CALL WriteDataToGroup(Hparams%sp%Jxy,paramsgrpId,"JXY")
		CALL WriteDataToGroup(Hparams%sp%h,paramsgrpId,"H")
		CALL WriteDataToGroup(Hparams%sp%gam,paramsgrpId,"GAMMA")
		CALL WriteDataToGroup(Hparams%sp%d,paramsgrpId,"D")
		CALL WriteDataToGroup(Hparams%sp%k,paramsgrpId,"K")
	CASE ('boson Hubbard')
		IF(.not.itp) THEN
			CALL WriteDataToGroup(Nmax,paramsgrpId,"Nmax")
		END IF
		CALL WriteDataToGroup(Hparams%bp%t,paramsgrpId,"T")
		CALL WriteDataToGroup(Hparams%bp%U,paramsgrpId,"U")
		CALL WriteDataToGroup(Hparams%bp%V,paramsgrpId,"V")
		CALL WriteDataToGroup(Hparams%bp%mu,paramsgrpId,"MU")
	CASE ('hardcore boson')
		CALL WriteDataToGroup(Hparams%hcbp%t,paramsgrpId,"T")
		CALL WriteDataToGroup(Hparams%hcbp%V,paramsgrpId,"V")
		CALL WriteDataToGroup(Hparams%hcbp%mu,paramsgrpId,"MU")
	CASE ('fermion Hubbard')
		CALL WriteDataToGroup(Hparams%fhp%t,paramsgrpId,"T")
		CALL WriteDataToGroup(Hparams%fhp%U,paramsgrpId,"U")
		CALL WriteDataToGroup(Hparams%fhp%V,paramsgrpId,"V")
		CALL WriteDataToGroup(Hparams%fhp%mu,paramsgrpId,"MU")
	CASE ('spinless fermions')
		CALL WriteDataToGroup(Hparams%sfp%t,paramsgrpId,"T")
		CALL WriteDataToGroup(Hparams%sfp%V,paramsgrpId,"V")
		CALL WriteDataToGroup(Hparams%sfp%mu,paramsgrpId,"MU")
	CASE DEFAULT
		PRINT *, "Hamiltonian type not recognized in WriteITPSettings!"
		PRINT *, "Use 'spin', 'boson Hubbard', 'hardcore boson', 'fermion Hubbard', or 'spinless fermions'."
		STOP 
END SELECT
CALL WriteDataToGroup(numQuenches,paramsgrpId,"numQuenches")
CALL WriteDataToGroup(rtpD%tau,paramsgrpId,"TAUS")
CALL WriteDataToGroup(rtpD%pow,paramsgrpId,"POWS")
CALL WriteDataToGroup(rtpD%gquench,paramsgrpId,"QUENCH IDS")
CALL WriteDataToGroup(rtpD%gi,paramsgrpId,"INITIAL QUENCH VALUES")
CALL WriteDataToGroup(rtpD%gf,paramsgrpId,"FINAL QUENCH VALUES")
CALL WriteDataToGroup(rtpD%nsteps,paramsgrpId,"NUM TIME STEPS")
CALL WriteDataToGroup(rtpD%stepsforStore,paramsgrpId,"STEPS FOR STORE")
CALL CloseGroup(paramsgrpID)
CALL CreateGroup("timesteps", timestepsgrpID, outputFileID) !Create timesteps group
CALL CloseGroup(timestepsgrpID)
CALL CloseHdf5File(outputFileID)


END SUBROUTINE WriteRTPSettings

SUBROUTINE WriteMeasures(Measures,HamiType,counter, time, truncErr, HparamsList, otherCounter)
!
!Purpose: Write the itp parameters into a group in the hdf5 file
!
IMPLICIT NONE
CHARACTER(len=*) :: HamiType
TYPE(Measure) :: Measures
INTEGER, INTENT(IN) :: counter, otherCounter
TYPE(HamiParamsList) :: HparamsList
INTEGER(HID_T) :: counterID,resultsID, energyID,entropiesID, avgID, localID, corrID, overlapID
INTEGER(HID_T) ::  qId, vnId,chainId, LEId, eID, timeID, teID, lpId
INTEGER(HID_T), ALLOCATABLE :: avIds(:), localIds(:), corrIds(:), hpId(:)
REAL(KIND=rKind), INTENT(IN) :: time, truncErr
CHARACTER(len=8) :: mString

CALL OpenHdf5File(outputName,outputFileID, openKind='rw')
CALL OpenGroup("timesteps", timestepsgrpID, outputFileID) !Open timesteps group

WRITE(mString,'(I8)') counter
CALL CreateSubGroup(mstring,counterID,timestepsgrpID) !Create new # subgroup
	CALL CreateSubGroup("results",resultsID,counterID) !Create new results subgroup
		CALL CreateSubGroup("cumulative truncation error",teID,resultsID) !Create new truncerr subgroup
			CALL WriteDataToMeanGroup(truncErr,teID,"value")
			CALL CloseGroup(teID)
		CALL CreateSubGroup("energy",energyID,resultsID) !Create new energy subgroup
			CALL WriteDataToMeanGroup(Measures%en,energyID,"value")
			CALL CloseGroup(energyID)
		CALL CreateSubGroup("Q measure",qId,resultsID) 
			CALL WriteDataToMeanGroup(REAL(Measures%ent%qme,KIND=rKind),qId,"value")
			CALL CloseGroup(qId)
		CALL CreateSubGroup("site entropy",vnId,resultsID) 
			CALL WriteDataToMeanGroup(REAL(Measures%ent%vN,KIND=rKind),vnId,"value")
			CALL CloseGroup(vnId)
		CALL CreateSubGroup("bond entropy",chainId,resultsID) 
			CALL WriteDataToMeanGroup(REAL(Measures%ent%chain,KIND=rKind),chainId,"value")
			CALL CloseGroup(chainId)
		CALL CreateSubGroup("Loschmidt Echo",LEId,resultsID) 
			CALL WriteDataToMeanGroup(ABS(Measures%overlap(1)%value)**2,LEId,"value")
			CALL CloseGroup(LEId)

SELECT CASE(HamiType)
	CASE ('spin')
		ALLOCATE(hpId(6))
		CALL CreateSubGroup("Local Props",lpId,counterID) !Create new props subgroup
			CALL WriteDataToGroup(time,lpId,"Time") !Write time data to /timesteps/#/
			CALL WriteDataToGroup(simId,lpId,"SIMID")
			CALL WriteDataToGroup(HparamsList%sp(Othercounter)%Jz,lpId,"Jz") 
			CALL WriteDataToGroup(HparamsList%sp(Othercounter)%Jxy,lpId,"Jxy") 
			CALL WriteDataToGroup(HparamsList%sp(Othercounter)%H,lpId,"H")
			CALL WriteDataToGroup(HparamsList%sp(Othercounter)%Gam,lpId,"Gamma") 
			CALL WriteDataToGroup(HparamsList%sp(Othercounter)%D,lpId,"D") 
			CALL WriteDataToGroup(HparamsList%sp(Othercounter)%K,lpId,"K") 
			CALL CloseGroup(lpId)
		CALL CreateSubGroup("Jz",hpId(1),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%sp(Othercounter)%Jz,hpId(1),"value")
			CALL CloseGroup(hpId(1))
		CALL CreateSubGroup("Jxy",hpId(2),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%sp(Othercounter)%Jxy,hpId(2),"value")
			CALL CloseGroup(hpId(2))
		CALL CreateSubGroup("H",hpId(3),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%sp(Othercounter)%H,hpId(3),"value")
			CALL CloseGroup(hpId(3))
		CALL CreateSubGroup("Gamma",hpId(4),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%sp(Othercounter)%Gam,hpId(4),"value")
			CALL CloseGroup(hpId(4))
		CALL CreateSubGroup("D",hpId(5),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%sp(Othercounter)%D,hpId(5),"value")
			CALL CloseGroup(hpId(5))
		CALL CreateSubGroup("K",hpId(6),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%sp(Othercounter)%K,hpId(6),"value")
			CALL CloseGroup(hpId(6))
		DEALLOCATE(hpId)

		ALLOCATE(avIds(4), localIds(4),corrIds(2))
		CALL CreateSubGroup("Avg Magnetization",avIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(1)%value,KIND=rKind),avIds(1),"value")
			CALL CloseGroup(avIds(1))
		CALL CreateSubGroup("Avg Magnetization^2",avIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(2)%value,KIND=rKind),avIds(2),"value")
			CALL CloseGroup(avIds(2))
		CALL CreateSubGroup("Avg Transverse Magnetization",avIds(3),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(3)%value,KIND=rKind),avIds(3),"value")
			CALL CloseGroup(avIds(3))
		CALL CreateSubGroup("Avg Transverse Magnetization^2",avIds(4),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(4)%value,KIND=rKind),avIds(4),"value")
			CALL CloseGroup(avIds(4))
		CALL CreateSubGroup("Local Magnetization",localIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(1)%value,KIND=rKind),localIds(1),"value")
			CALL CloseGroup(localIds(1))
		CALL CreateSubGroup("Local Magnetization^2",localIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(2)%value,KIND=rKind),localIds(2),"value")
			CALL CloseGroup(localIds(2))
		CALL CreateSubGroup("Local Transverse Magnetization",localIds(3),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(3)%value,KIND=rKind),localIds(3),"value")
			CALL CloseGroup(localIds(3))
		CALL CreateSubGroup("Local Transverse Magnetization^2",localIds(4),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(4)%value,KIND=rKind),localIds(4),"value")
			CALL CloseGroup(localIds(4))
		CALL CreateSubGroup("<Sz Sz>",corrIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%corr(1)%value,KIND=rKind),corrIds(1),"value")
			CALL CloseGroup(corrIds(1))
		CALL CreateSubGroup("<Sx Sx>",corrIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%corr(2)%value,KIND=rKind),corrIds(2),"value")
			CALL CloseGroup(corrIds(2))
		DEALLOCATE(avIds, localIds, corrIds)
	CASE ('boson Hubbard')
		ALLOCATE(hpId(4))
		CALL CreateSubGroup("Local Props",lpId,counterID) !Create new props subgroup
			CALL WriteDataToGroup(time,lpId,"Time") !Write time data to /timesteps/#/
			CALL WriteDataToGroup(simId,lpId,"SIMID")
			CALL WriteDataToGroup(HparamsList%bp(Othercounter)%mu,lpId,"mu") 
			CALL WriteDataToGroup(HparamsList%bp(Othercounter)%V,lpId,"V") 
			CALL WriteDataToGroup(HparamsList%bp(Othercounter)%t,lpId,"t")
			CALL WriteDataToGroup(HparamsList%bp(Othercounter)%U,lpId,"U") 
			CALL CloseGroup(lpId)
		CALL CreateSubGroup("mu",hpId(1),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%bp(Othercounter)%mu,hpId(1),"value")
			CALL CloseGroup(hpId(1))
		CALL CreateSubGroup("V",hpId(2),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%bp(Othercounter)%V,hpId(2),"value")
			CALL CloseGroup(hpId(2))
		CALL CreateSubGroup("t",hpId(3),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%bp(Othercounter)%t,hpId(3),"value")
			CALL CloseGroup(hpId(3))
		CALL CreateSubGroup("U",hpId(4),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%bp(Othercounter)%U,hpId(4),"value")
			CALL CloseGroup(hpId(4))
		DEALLOCATE(hpId)
		ALLOCATE(avIds(2), localIds(2),corrIds(2))
		CALL CreateSubGroup("Avg Number",avIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(1)%value,KIND=rKind),avIds(1),"value")
			CALL CloseGroup(avIds(1))
		CALL CreateSubGroup("Avg Number^2",avIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(2)%value,KIND=rKind),avIds(2),"value")
			CALL CloseGroup(avIds(2))
		CALL CreateSubGroup("Local Number",localIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(1)%value,KIND=rKind),localIds(1),"value")
			CALL CloseGroup(localIds(1))
		CALL CreateSubGroup("Local Number^2",localIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(2)%value,KIND=rKind),localIds(2),"value")
			CALL CloseGroup(localIds(2))
		CALL CreateSubGroup("<n n>",corrIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%corr(1)%value,KIND=rKind),corrIds(1),"value")
			CALL CloseGroup(corrIds(1))
		CALL CreateSubGroup("<adag a>",corrIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%corr(2)%value,KIND=rKind),corrIds(2),"value")
			CALL CloseGroup(corrIds(2))
		DEALLOCATE(avIds, localIds, corrIds)
	CASE ('hardcore boson')
		ALLOCATE(hpId(3))
		CALL CreateSubGroup("Local Props",lpId,counterID) !Create new props subgroup
			CALL WriteDataToGroup(time,lpId,"Time") !Write time data to /timesteps/#/
			CALL WriteDataToGroup(simId,lpId,"SIMID")
			CALL WriteDataToGroup(HparamsList%hcbp(Othercounter)%mu,lpId,"mu") 
			CALL WriteDataToGroup(HparamsList%hcbp(Othercounter)%V,lpId,"V") 
			CALL WriteDataToGroup(HparamsList%hcbp(Othercounter)%t,lpId,"t")
			CALL CloseGroup(lpId)
		CALL CreateSubGroup("mu",hpId(1),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%hcbp(Othercounter)%mu,hpId(1),"value")
			CALL CloseGroup(hpId(1))
		CALL CreateSubGroup("V",hpId(2),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%hcbp(Othercounter)%V,hpId(2),"value")
			CALL CloseGroup(hpId(2))
		CALL CreateSubGroup("t",hpId(3),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%hcbp(Othercounter)%t,hpId(3),"value")
			CALL CloseGroup(hpId(3))
		DEALLOCATE(hpId)
		ALLOCATE(avIds(2), localIds(2),corrIds(2))
		CALL CreateSubGroup("Avg Number",avIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(1)%value,KIND=rKind),avIds(1),"value")
			CALL CloseGroup(avIds(1))
		CALL CreateSubGroup("Avg Number^2",avIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(2)%value,KIND=rKind),avIds(2),"value")
			CALL CloseGroup(avIds(2))
		CALL CreateSubGroup("Local Number",localIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(1)%value,KIND=rKind),localIds(1),"value")
			CALL CloseGroup(localIds(1))
		CALL CreateSubGroup("Local Number^2",localIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(2)%value,KIND=rKind),localIds(2),"value")
			CALL CloseGroup(localIds(2))
		CALL CreateSubGroup("<n n>",corrIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%corr(1)%value,KIND=rKind),corrIds(1),"value")
			CALL CloseGroup(corrIds(1))
		CALL CreateSubGroup("<adag a>",corrIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%corr(2)%value,KIND=rKind),corrIds(2),"value")
			CALL CloseGroup(corrIds(2))
		DEALLOCATE(avIds, localIds, corrIds)
	CASE ('fermion Hubbard')
		ALLOCATE(hpId(4))
		CALL CreateSubGroup("Local Props",lpId,counterID) !Create new props subgroup
			CALL WriteDataToGroup(time,lpId,"Time") !Write time data to /timesteps/#/
			CALL WriteDataToGroup(simId,lpId,"SIMID")
			CALL WriteDataToGroup(HparamsList%fhp(Othercounter)%mu,lpId,"mu") 
			CALL WriteDataToGroup(HparamsList%fhp(Othercounter)%V,lpId,"V") 
			CALL WriteDataToGroup(HparamsList%fhp(Othercounter)%t,lpId,"t")
			CALL WriteDataToGroup(HparamsList%fhp(Othercounter)%U,lpId,"U") 
			CALL CloseGroup(lpId)
		CALL CreateSubGroup("mu",hpId(1),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%fhp(Othercounter)%mu,hpId(1),"value")
			CALL CloseGroup(hpId(1))
		CALL CreateSubGroup("V",hpId(2),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%fhp(Othercounter)%V,hpId(2),"value")
			CALL CloseGroup(hpId(2))
		CALL CreateSubGroup("t",hpId(3),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%fhp(Othercounter)%t,hpId(3),"value")
			CALL CloseGroup(hpId(3))
		CALL CreateSubGroup("U",hpId(4),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%fhp(Othercounter)%U,hpId(4),"value")
			CALL CloseGroup(hpId(4))
		DEALLOCATE(hpId)
		ALLOCATE(avIds(4), localIds(4),corrIds(5))
		CALL CreateSubGroup("Avg Number",avIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(1)%value,KIND=rKind),avIds(1),"value")
			CALL CloseGroup(avIds(1))
		CALL CreateSubGroup("Avg Number^2",avIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(2)%value,KIND=rKind),avIds(2),"value")
			CALL CloseGroup(avIds(2))
		CALL CreateSubGroup("Avg Magnetization",avIds(3),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(3)%value,KIND=rKind),avIds(3),"value")
			CALL CloseGroup(avIds(3))
		CALL CreateSubGroup("Avg Magnetization^2",avIds(4),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(4)%value,KIND=rKind),avIds(4),"value")
			CALL CloseGroup(avIds(4))
		CALL CreateSubGroup("Local Number",localIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(1)%value,KIND=rKind),localIds(1),"value")
			CALL CloseGroup(localIds(1))
		CALL CreateSubGroup("Local Number^2",localIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(2)%value,KIND=rKind),localIds(2),"value")
			CALL CloseGroup(localIds(2))
		CALL CreateSubGroup("Local Magnetization",localIds(3),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(3)%value,KIND=rKind),localIds(3),"value")
			CALL CloseGroup(localIds(3))
		CALL CreateSubGroup("Local Magnetization^2",localIds(4),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(4)%value,KIND=rKind),localIds(4),"value")
			CALL CloseGroup(localIds(4))
		CALL CreateSubGroup("<n n>",corrIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%corr(1)%value,KIND=rKind),corrIds(1),"value")
			CALL CloseGroup(corrIds(1))
		CALL CreateSubGroup("<Sz Sz>",corrIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%corr(2)%value,KIND=rKind),corrIds(2),"value")
			CALL CloseGroup(corrIds(2))
!		CALL CreateSubGroup("Pairing",corrIds(3),resultsID)
!			CALL WriteDataToMeanGroup(REAL(Measures%corr(3)%value,KIND=rKind),corrIds(3),"value")
!			CALL CloseGroup(corrIds(3))
		CALL CreateSubGroup("<adag_up a_up>",corrIds(4),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%fermicorr(1)%value,KIND=rKind),corrIds(4),"value")
			CALL CloseGroup(corrIds(4))
		CALL CreateSubGroup("<adag_down a_down>",corrIds(5),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%fermicorr(2)%value,KIND=rKind),corrIds(5),"value")
			CALL CloseGroup(corrIds(5))
		DEALLOCATE(avIds, localIds, corrIds)
	CASE ('spinless fermions')
		ALLOCATE(hpId(3))
		CALL CreateSubGroup("Local Props",lpId,counterID) !Create new props subgroup
			CALL WriteDataToGroup(time,lpId,"Time") !Write time data to /timesteps/#/
			CALL WriteDataToGroup(simId,lpId,"SIMID")
			CALL WriteDataToGroup(HparamsList%sfp(Othercounter)%mu,lpId,"mu") 
			CALL WriteDataToGroup(HparamsList%sfp(Othercounter)%V,lpId,"V") 
			CALL WriteDataToGroup(HparamsList%sfp(Othercounter)%t,lpId,"t")
			CALL CloseGroup(lpId)
		CALL CreateSubGroup("mu",hpId(1),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%sfp(Othercounter)%mu,hpId(1),"value")
			CALL CloseGroup(hpId(1))
		CALL CreateSubGroup("V",hpId(2),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%sfp(Othercounter)%V,hpId(2),"value")
			CALL CloseGroup(hpId(2))
		CALL CreateSubGroup("t",hpId(3),resultsID)
			CALL WriteDataToMeanGroup(HparamsList%sfp(Othercounter)%t,hpId(3),"value")
			CALL CloseGroup(hpId(3))
		DEALLOCATE(hpId)

		ALLOCATE(avIds(2), localIds(2),corrIds(2))
		CALL CreateSubGroup("Avg Number",avIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(1)%value,KIND=rKind),avIds(1),"value")
			CALL CloseGroup(avIds(1))
		CALL CreateSubGroup("Avg Number^2",avIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%avg(2)%value,KIND=rKind),avIds(2),"value")
			CALL CloseGroup(avIds(2))
		CALL CreateSubGroup("Local Number",localIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(1)%value,KIND=rKind),localIds(1),"value")
			CALL CloseGroup(localIds(1))
		CALL CreateSubGroup("Local Number^2",localIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%local(2)%value,KIND=rKind),localIds(2),"value")
			CALL CloseGroup(localIds(2))
		CALL CreateSubGroup("<n n>",corrIds(1),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%corr(1)%value,KIND=rKind),corrIds(1),"value")
			CALL CloseGroup(corrIds(1))
		CALL CreateSubGroup("<adag a>",corrIds(2),resultsID)
			CALL WriteDataToMeanGroup(REAL(Measures%fermicorr(1)%value,KIND=rKind),corrIds(2),"value")
			CALL CloseGroup(corrIds(2))
		DEALLOCATE(avIds, localIds, corrIds)
	CASE DEFAULT
		PRINT *, "Hamiltonian type not recognized in WriteITPSettings!"
		PRINT *, "Use 'spin', 'boson Hubbard', 'hardcore boson', 'fermion Hubbard', or 'spinless fermions'."
		STOP 
END SELECT
CALL CloseGroup(resultsId)
CALL CloseGroup(counterID)
CALL CloseGroup(timeStepsgrpID)
CALL CloseHdf5File(outputFileID)


END SUBROUTINE WriteMeasures





END MODULE Hdf5Interface
