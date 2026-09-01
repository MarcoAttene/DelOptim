ECHO OFF

WHERE /Q cmake
IF %ERRORLEVEL% NEQ 0 (
	ECHO ERROR! Could not find CMake in your system path. Please verify your CMake installation.
	EXIT /B 1
)

mkdir build
cd build
cmake ..

IF %ERRORLEVEL% NEQ 0 (
	ECHO ERROR! CMake failed.
	EXIT /B 1
)

cmake --build . --config Release

IF %ERRORLEVEL% NEQ 0 (
	ECHO ERROR! Building failed.
	EXIT /B 1
)

cd Release
delmesher -huwxyz ../../input_models/boeing_part.off

echo Files chamfered_plc.off , DR_interface.off , enrichedCDT_constrainedFaces.off , DR_mesh.tet and enrichedCDT_mesh.tet have been created in directory build/Release

set /p DUMMY=Hit ENTER to terminate this script...