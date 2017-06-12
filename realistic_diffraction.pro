TARGET = diffraction_huygens
TEMPLATE = app

#QMAKE_LFLAGS += -framework OpenCL
#QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.9

SOURCES += main.cpp \
    huygens.cpp \
    PFMReadWrite.cpp \
    integration.cpp


HEADERS  += huygens.h \
    PFMReadWrite.h \
    integration.h

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
