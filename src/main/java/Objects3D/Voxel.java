/**
 *
 *  Coordinates3D v1, 1 nov. 2009
    Fabrice P Cordelieres, fabrice.cordelieres at gmail.com

    Copyright (C) 2009 Fabrice P. Cordelieres

    License:
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

package Objects3D;

public class Voxel extends Coordinates3D{
    public int val;
    public boolean isSurf;

    public Voxel(int x, int y, int z, int val, boolean isSurf){
        this.x=x;
        this.y=y;
        this.z=z;
        this.val=val;
        this.isSurf=isSurf;
    }
    
    public Voxel(Coordinates3D coord3D, int val, boolean isSurf){
        x=coord3D.x;
        y=coord3D.y;
        z=coord3D.z;
        this.val=val;
        this.isSurf=isSurf;
    }



    public Coordinates3D getCoodinates3D(){
        return new Coordinates3D(x, y, z);
    }

}
