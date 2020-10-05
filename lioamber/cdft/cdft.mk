INCLUDES := cdft_data.o
OBJECTS += $(INCLUDES)

$(OBJPATH)/cdft_subs.o: cdft_data.o cdft_main.f90 cdft_input.f90 \
                        cdft_init_fin.f90 cdft_utils.f90         \
								cdft_mixed_utils.f90
$(OBJPATH)/cdft_subs.o: $(INCLUDES) cdft.mk
