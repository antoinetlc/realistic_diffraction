QT       += core gui
QT       += opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = diffraction_huygens
TEMPLATE = app

#QMAKE_LFLAGS += -framework OpenCL
#QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.9

SOURCES += main.cpp \
    huygens.cpp \
    mathfunctions.cpp \
    PFMReadWrite.cpp


HEADERS  += huygens.h \
    mathfunctions.h \
    PFMReadWrite.h

FORMS    +=
OTHER_FILES +=

DISTFILES +=


INCLUDEPATH += "C:/glew-1.13.0/include"
LIBS += "C:/glew-1.13.0/lib/Release/x64/glew32.lib"

##################### CUDA   ##############################


win32:{

    #.cu files compiled with NVCC
    #OBJECTS += cuda/vectoradd_cuda.obj

    #Library directories include
    INCLUDEPATH += "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.0\include"

    #Library directories
    QMAKE_LIBDIR += "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.0\lib\x64"

    LIBS +=  "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.0\lib\x64\cuda.lib"
    LIBS +=  "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.0\lib\x64\cudart.lib"
    LIBS +=  "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.0\lib\x64\cufft.lib"

}
else:unix{
    SYSTEM_NAME = unix
    CUDA_DIR = /usr/local/cuda-7.0/
    INCLUDEPATH += $$CUDA_DIR/include
    QMAKE_LIBDIR += $$CUDA_DIR/lib64
    LIBS += -lcudart -lcuda -lcufft
    CUDA_ARCH     = sm_11 #Compute capability 3.5 on GTX 630

    # NVCC flags (set by default)
    NVCCFLAGS     = --compiler-options -fno-strict-aliasing -use_fast_math --ptxas-options=-v

    # Prepare the extra compiler configuration (taken from the nvidia forum)
    CUDA_INC = $$join(INCLUDEPATH,' -I','-I',' ')

    cuda.commands = $$CUDA_DIR/bin/nvcc -m64 -O3 -arch=$$CUDA_ARCH -c $$NVCCFLAGS \
                    $$CUDA_INC $$LIBS  ${QMAKE_FILE_NAME} -o ${QMAKE_FILE_OUT} \
                    2>&1 | sed -r \"s/\\(([0-9]+)\\)/:\\1/g\" 1>&2
    # nvcc error printout format ever so slightly different from gcc
    # http://forums.nvidia.com/index.php?showtopic=171651

    cuda.dependency_type = TYPE_C
    cuda.depend_command = $$CUDA_DIR/bin/nvcc -O3 -M $$CUDA_INC $$NVCCFLAGS   ${QMAKE_FILE_NAME}

    cuda.input = CUDA_SOURCES
    cuda.output = ${OBJECTS_DIR}${QMAKE_FILE_BASE}_cuda.o

    # Qt will add cuda to the Makefile
    QMAKE_EXTRA_COMPILERS += cuda

}

##################### OpenCV   ##############################



win32:{
    CONFIG(debug, debug|release)
    {

         INCLUDEPATH += "C:\\OpenCV2411\\build\\include"

         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_core2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_highgui2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_imgproc2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_features2d2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_calib3d2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_contrib2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_flann2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_gpu2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_legacy2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_ml2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_nonfree2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_objdetect2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_ocl2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_photo2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_stitching2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_superres2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_ts2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_video2411.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_videostab2411.lib"

    }
    CONFIG(release, debug|release)
    {
         INCLUDEPATH += "C:\\OpenCV2411\\build\\include"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_core2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_highgui2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_imgproc2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_features2d2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_calib3d2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_contrib2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_flann2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_gpu2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_legacy2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_ml2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_nonfree2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_objdetect2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_ocl2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_photo2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_stitching2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_superres2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_ts2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_video2411d.lib"
         LIBS += "C:\\OpenCV2411\\build\\x64\\vc12\\lib\\opencv_videostab2411d.lib"
    }
}
else:unix{
    INCLUDEPATH += /usr/local/include/
    LIBS += -L/usr/local/lib
    LIBS += -lopencv_core
    LIBS += -lopencv_imgproc
    LIBS += -lopencv_highgui
    LIBS += -lopencv_ml
    LIBS += -lopencv_video
    LIBS += -lopencv_features2d
    LIBS += -lopencv_calib3d
    LIBS += -lopencv_objdetect
    LIBS += -lopencv_contrib
    LIBS += -lopencv_legacy
    LIBS += -lopencv_flann
}

##################### DLib   ##############################

win32:{
INCLUDEPATH += C:\Libraries\dlib-18.16
}
else:unix{

}

##################### OpenCV  WINDOWS MINGW 32 bits ##############################

#INCLUDEPATH += C:\OpenCV-mingw\install\include
#LIBS += -LC:\OpenCV-mingw\install\x64\mingw\bin \
#    libopencv_core2411 \
#    libopencv_highgui2411 \
#    libopencv_imgproc2411 \
#    libopencv_features2d2411 \
#    libopencv_calib3d2411
