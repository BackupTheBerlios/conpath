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

include ../Make.inc

install: 
	$(MAKE) -f Makefile.preprocess source_list
	$(MAKE) -f Makefile.preprocess all
	$(MAKE) -f Makefile.comp config
	$(MAKE) -f Makefile.comp install

clean:
	$(MAKE) -f Makefile.preprocess clean
	$(MAKE) -f Makefile.comp clean

allclean:
	$(MAKE) -f Makefile.preprocess allclean
	$(MAKE) -f Makefile.comp allclean

config:
	$(MAKE) -f Makefile.preprocess source_list
	$(MAKE) -f Makefile.preprocess all
	$(MAKE) -f Makefile.comp config

purge:
	 @find . -name "*~" -exec rm -f  {} \;

# End of ConPath Makefile ...
