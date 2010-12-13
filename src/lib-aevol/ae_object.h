//*****************************************************************************
// S.P.E.A.R. library - Simulator of Physical Environment for Animat Research
// Copyright (C) 2003  S�bastien GRIPON, Fran�ois PERRIN, H�di SOULA - PRISMa
// Web: http://prisma.insa-lyon.fr
// Original Authors : Guillaume BESLON, H�di SOULA
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//*****************************************************************************

// File Information
// $Id: f_object.h,v 1.3 2005/05/03 13:11:28 cknibbe Exp $

#ifndef __AE_OBJECT_H__
#define __AE_OBJECT_H__

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>


class ae_object
{
  public :
    ae_object( void ){};
    virtual ~ae_object( void ){};
    
  protected :
    ae_object( const ae_object &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

};

#endif // __AE_OBJECT_H__
