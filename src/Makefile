###############################################
#       Makefile for Quality Assessment       #
###############################################

#### CFLAGS     	Flags to pass to C compiler.
CFLAGS          	= -O2 -lm -Wall

all:			calc_neff  calculate_TA_prediction  convert  gentable parse_metapsicov parse_proxy_map

clean:
			rm -rf *.o calc_neff  calculate_TA_prediction  convert  gentable parse_metapsicov parse_proxy_map

install: 
			cp calc_neff calculate_TA_prediction convert gentable parse_metapsicov parse_proxy_map ../bin


calc_neff:			calc_neff.c
	$(CC) $(CFLAGS) calc_neff.c -o calc_neff

calculate_TA_prediction:          calculate_TA_prediction.c
	$(CC) $(CFLAGS) calculate_TA_prediction.c -o calculate_TA_predicion

convert:          Converter.c
	$(CC) $(CFLAGS) Converter.c -o convert

gentable:          generate_table.c
	$(CC) $(CFLAGS) generate_table.c -o gentable

parse_metapsicov:          parse_metapsicov.c
	$(CC) $(CFLAGS) parse_metapsicov.c -o parse_metapsicov

parse_proxy_map:          parse_proxy_map.c
	$(CC) $(CFLAGS) parse_proxy_map.c -o parse_proxy_map


