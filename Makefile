USER = $(shell whoami)
MAKE = /usr/bin/make
DATA_DIR = data
SRC = src
TEX_DIR = /Users/$(USER)/Dropbox/fluxnet_heatwaves_paper/figures
FIG_DIR = /Users/$(USER)/Dropbox/fluxnet_heatwaves_paper/figures/figs
SUB_FIG_DIR = /Users/$(USER)/Dropbox/fluxnet_heatwaves_paper/submission_figures

all: data figs paper
data: $(DATA_DIR)/fluxnet2015_all_events.csv \
	  $(DATA_DIR)/ozflux_all_events.csv \
	  $(DATA_DIR)/ozflux_slopes_GPP.csv \
	  $(DATA_DIR)/ozflux_slopes_Qle.csv \
	  $(DATA_DIR)/fluxnet2015_heatwaves_gpp_sqrt_d_et.csv \
	  $(DATA_DIR)/fluxnet2015_heatwaves_gpp_sqrt_d_et_non.csv
figs: $(FIG_DIR)/all_events_GPP_positive.pdf \
		$(FIG_DIR)/all_events_GPP_negative.pdf \
		$(FIG_DIR)/all_events_Qle_positive.pdf \
		$(FIG_DIR)/all_events_Qle_negative.pdf \
		$(FIG_DIR)/ozflux_histogram_GPP_ozflux.pdf \
		$(FIG_DIR)/ozflux_histogram_Qle_ozflux.pdf \
		$(FIG_DIR)/ozflux_heatwave_vs_non_heatwave_gppsqrtd_et.pdf \
		$(FIG_DIR)/timing_Qle.pdf \
		$(FIG_DIR)/all_events_GPP_FLUXNET_negative.pdf \
		$(FIG_DIR)/all_events_Qle_FLUXNET_positive.pdf
paper: $(TEX_DIR)/figures.pdf

##
# Make data files
##
$(DATA_DIR)/fluxnet2015_all_events.csv:	$(SRC)/get_fluxes_for_all_events_flux_data.py
	python $<
$(DATA_DIR)/ozflux_all_events.csv:	$(SRC)/get_fluxes_for_all_events_flux_data.py
	python $<
$(DATA_DIR)/ozflux_heatwaves_gpp_sqrt_d_et_non.csv:	$(SRC)/get_gpp_sqrt_D_et_for_heatwaves_flux_data.py
	python $<
$(DATA_DIR)/ozflux_heatwaves_gpp_sqrt_d_et.csv:	$(SRC)/get_gpp_sqrt_D_et_for_heatwaves_flux_data.py
	python $<

##
# Figures - set these up so that if any data changes, they are *ALL* remade
##

$(FIG_DIR)/all_events_GPP_negative.pdf:	$(SRC)/plot_GPP_at_all_events_above_Tthreh.py \
										$(DATA_DIR)/ozflux_all_events.csv
	python $<
$(FIG_DIR)/all_events_GPP_positive.pdf:	$(SRC)/plot_GPP_at_all_events_above_Tthreh.py \
										$(DATA_DIR)/ozflux_all_events.csv
	python $<
$(FIG_DIR)/all_events_Qle_positive.pdf:	$(SRC)/plot_Qle_at_all_events_above_Tthreh.py \
										$(DATA_DIR)/ozflux_all_events.csv
	python $<
$(FIG_DIR)/all_events_Qle_negative.pdf:	$(SRC)/plot_Qle_at_all_events_above_Tthreh.py \
										$(DATA_DIR)/ozflux_all_events.csv
	python $<
$(FIG_DIR)/ozflux_histogram_GPP_ozflux.pdf:	$(SRC)/plot_ozflux_slope_histogram_GPP.py \
											$(DATA_DIR)/ozflux_slopes_GPP.csv
	python $<
$(FIG_DIR)/ozflux_histogram_Qle_ozflux.pdf:	$(SRC)/plot_ozflux_slope_histogram_Qle.py \
											$(DATA_DIR)/ozflux_slopes_Qle.csv
	python $<
$(FIG_DIR)/ozflux_heatwave_vs_non_heatwave_gppsqrtd_et.pdf:	$(SRC)/plot_ozflux_gpp_sqrtd_et.py \
															$(DATA_DIR)/ozflux_heatwaves_gpp_sqrt_d_et_non.csv \
															$(DATA_DIR)/ozflux_heatwaves_gpp_sqrt_d_et.csv \
															$(DATA_DIR)/fluxnet2015_heatwaves_gpp_sqrt_d_et_non.csv \
															$(DATA_DIR)/fluxnet2015_heatwaves_gpp_sqrt_d_et.csv
	python $<
$(FIG_DIR)/timing_Qle.pdf:	$(SRC)/plot_timing_of_LE_slopes.py \
							$(DATA_DIR)/ozflux_slopes_Qle.csv
	python $<
$(FIG_DIR)/all_events_GPP_FLUXNET_negative.pdf:	$(SRC)/plot_GPP_at_all_events_above_Tthreh_FLUXNET.py \
												$(DATA_DIR)/ozflux_all_events.csv
	python $<
$(FIG_DIR)/all_events_Qle_FLUXNET_positive.pdf:	$(SRC)/plot_Qle_at_all_events_above_Tthreh_FLUXNET.py \
												$(DATA_DIR)/ozflux_all_events.csv
	python $<


$(TEX_DIR)/figures.pdf:	$(TEX_DIR)/figures.tex \
						$(FIG_DIR)/all_events_GPP_negative.pdf \
						$(FIG_DIR)/all_events_GPP_positive.pdf \
						$(FIG_DIR)/all_events_Qle_positive.pdf \
						$(FIG_DIR)/all_events_Qle_negative.pdf \
						$(FIG_DIR)/ozflux_histogram_GPP_ozflux.pdf \
						$(FIG_DIR)/ozflux_histogram_Qle_ozflux.pdf \
						$(FIG_DIR)/ozflux_heatwave_vs_non_heatwave_gppsqrtd_et.pdf \
						$(FIG_DIR)/timing_Qle.pdf \
						$(FIG_DIR)/all_events_GPP_FLUXNET_negative.pdf \
						$(FIG_DIR)/all_events_Qle_FLUXNET_positive.pdf
	cd $(TEX_DIR) && $(MAKE) clean && $(MAKE)
	cd $(SUB_FIG_DIR) && ./get_figs.BASH && $(MAKE) clean && $(MAKE)
cleanfigs:
	rm -f $(FIG_DIR)/*.pdf

.PHONY : clean

clean:
	rm -f $(DATA_DIR)/*.csv $(FIG_DIR)/*.pdf
