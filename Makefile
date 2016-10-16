PLUGIN_NAME = phase_response_curve

RTXI_INCLUDES =

HEADERS = phase-response-curve.h\

SOURCES = phase-response-curve.cpp \
          moc_phase-response-curve.cpp\

LIBS = -lqwt -lrtplot

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
