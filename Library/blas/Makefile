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
#     Makefile to compile BLAS routines

include ../../Make.inc
include objlist

blas: $(OBJS)
	
clean:
	rm -f $(OBJS) 

config: 
	$(CONPATH_BIN)/Config

.SUFFIXES: $(SUFFIXES) .f

.f.o:
	$(F90) $(F90LAPACKFLAGS) -c $<
#
# All dependencies beetwen modules and files have been stored in 
# dependencies file
#

include dependencies

# End of ConPath Makefile ...
