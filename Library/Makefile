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
#     Makefile to compile Library needed by Conpath Program

include ../Make.inc

library: 
	@cd lapack; $(MAKE) lapack
	@cd blas; $(MAKE) blas  
	rm -f $(LAPACK30)
	@set OBJ=`find . -name "*.o" -exec echo {} \;`
	ar -cru $(LAPACK30) ./*/*o ./*/*/*o
	ranlib $(LAPACK30)

clean:
	@cd lapack; $(MAKE) clean
	@cd blas; $(MAKE) clean
	rm -f $(LAPACK30)

config:
	@cd lapack; $(MAKE) config
	@cd blas; $(MAKE) config

# End of ConPath Makefile ...
