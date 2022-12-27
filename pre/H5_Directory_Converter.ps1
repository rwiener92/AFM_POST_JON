# H5_Directory_Converter.ps1
# This script will convert all ARDF files in directory to h5 files with the same name
# Must first be in directory which has your ARDF files
### NOTE: do in powershell script, since WSL2 bash script wont run .exe


# ARDFtoHDF.exe location = C:\Users\Rob\Documents\GitHub\AFM_POST_JON\ARDFtoHDF5.exe
# ARDFtoHDF.exe <A.ARDF> <B.h5>

$dirfls = $(ls -Filter *.ARDF)
# ls | Where-Object Name .ARDF -Match

Write-Host ($dirfls).count items to convert
#$dirfls | Measure-Object


foreach ($item in $dirfls){

$newname = $item -replace "\.ARDF", ".h5"
Write-Host $item $newname
C:\Users\Rob\Documents\GitHub\AFM_POST_JON\ARDFtoHDF5.exe $item $newname

}
Write-Host End of foreach loop.



