#
# Program ConPath - Generical program to perform dynamics on a single state
#                 - Perform dynamics into a conical intersection
#                 - Perform dynamics in internal coordinates
#
#
#     Copyright (C) 2002-2004  Teodoro Laino                                
#                              Daniele Passerone
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
#     Makefile to compile all ConPath
#

include Make.inc

all: 
	@cd Library; $(MAKE) library
	@cd src; $(MAKE) install

clean:
	@cd Library; $(MAKE) clean
	@cd src; $(MAKE) allclean
	@find . -name "work.pc*" -exec rm -f  {} \;	

config:
	@cd Library; $(MAKE) config
	@cd src; $(MAKE) config

purge: 
	@cd Library; $(MAKE) purge
	@cd src; $(MAKE) purge

# End of ConPath Makefile ...
