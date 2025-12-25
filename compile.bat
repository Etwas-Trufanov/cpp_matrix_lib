@echo off
meson setup build
cd build
ninja
copy cpplibtest.exe ..\
cd ..
rmdir /s /q build
echo Done!
