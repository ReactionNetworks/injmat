/* Copyright (C) 2013, Murad Banaji
 *
 * This file is part of QUALMAT
 *
 * QUALMAT is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2, 
 * or (at your option) any later version.
 *
 * QUALMAT is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUALMAT: see the file COPYING.  If not, write to 
 * the Free Software Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA. 

 * For comments, queries, or bug-reports about QUALMAT, please 
 * email qualmat@gmail.com

 */

#include "crn.h"

int main(int argc, char *argv[])
{
  int q=1; // q=1 is quick, q=0 is slow but more info
  if(argc<2){
    fprintf(stderr, "No data file supplied, using default file \"datfiles/defreac\".\n\n");
    analysereacs("datfiles/defreac", q);
  }
  else
    analysereacs(argv[1], q); // quick

  return 0;

}

