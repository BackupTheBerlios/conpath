# Makefile include file .. Platform dependent
#
# LINUX-IFC-OS
#
# Program ConPath - Generical program to perform dynamics 
#                   in search of Conical Intersections (and not only!..)
#
#     Copyright (C) 2002  Teodoro Laino                                
#                         Daniele Passerone
#                                                                      
#     This program is free software; you can redistribute it and/or    
#     modify it under the terms of the GNU General Public License      
#     as published by the Free Software Foundation; either version 2   
#     of the License, or (at your option) any later version.          
#     This program is distributed in the hope that it will be useful, 
#     but WITHOUT ANY WARRANTY; without even the implied warranty of  
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   
#     GNU General Public License for more details.                    
#                                                                     
#     For information please contact author at: t.laino@sns.it        

CONPATH_BIN=$(CONPATH_DIR)/bin
MODDIR=$(CONPATH_DIR)/src/modules/
CPP = /usr/bin/cpp
CPPFLAGS = -O -P -traditional -DIFC
CC = cc
CFLAGS = -O
FC = /usr/sepp/bin/ifc-7.1-ms
FFLAGS = -O3  -Vaxlib -tpp7
F90 = /usr/sepp/bin/ifc-7.1-ms
F90FLAGS = -O3  -Vaxlib -tpp7  \
	    -module $(MODDIR)  -I$(MODDIR) 
F90LAPACKFLAGS = -O3 -Vaxlib  -tpp7 -w \
	    -module $(MODDIR)  -I$(MODDIR) 
LDFLAGS = $(F90FLAGS) 
LIBFLAGS= -L$(CONPATH_DIR)/lib -l$(LIBS)

# End of Makefile include file - Conpath


# Include Special for Lapack3.0

OPTS3    = -c $(F90FLAGS)  
OPTL     = -o
LAPACK30 = $(CONPATH_DIR)/lib/liblapack30.a

LIBS     = lapack30
SUF      = f90

# End of Makefile include Lapack30
