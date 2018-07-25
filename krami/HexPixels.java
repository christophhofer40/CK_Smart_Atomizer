package krami;
import ij.*; //needed for ij.Macro

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;

public class HexPixels
{
	private int A = -1;
	private int W = -1;
	private int R = -1;
	private int cq = 0;
	private int cr = 0;
	private int cx = cq - (cr - (cr&1)) / 2;
	private int cz = cr;
	private int cy = -cx - cz;
	private boolean[] hex_mask = null;
	
	public int hpMax = -1;
	public int[] hp = null;
	
	public HexPixels(int r)
	{
		R = r;
		W = 2*r;
		A = W*W;
		cq = r;
		cr = r;
		cx = cq - (cr - (cr&1)) / 2;
		cz = cr;
		cy = -cx - cz;
		hex_mask = new boolean[A];
		hpMax = ( 3 * A) / 4; // pixels per hex domain
		hp = new int[hpMax];
		int hpC = 0;
		for(int qr = 0; qr < A; ++qr)
		{
			final boolean ins = inside(qr);
			if(ins)
			{
				hex_mask[qr] = true;
				hp[hpC++] = qr;
			}
		}
		//System.out.println("created new HexPixels("+ r +")");
	}

	public int[] to_xyz(int qr)
	{
		final int q = qr % W;
		final int r = qr / W;
		final int x = q - (r - (r&1)) / 2 - cx;
		final int z = r - cz;
		final int y = -x - z;
		return new int[]{x,y,z};
	}
	
	public int get_Width()
	{	return W;}
	
	public int radius(int qr)
	{
		final int q = qr % W;
		final int r = qr / W;
		final int x = q - (r - (r&1)) / 2 - cx;
		final int z = r - cz;
		final int y = -x - z;
		return (Math.abs(x) + Math.abs(y) + Math.abs(z))/2;
	}
	
	public int hexdist(int q1, int r1, int q2, int r2)
	{
		final int x1 = q1 - (r1 - (r1&1)) / 2;
		final int x2 = q2 - (r2 - (r2&1)) / 2;
		final int dx = x2 - x1;
		final int dz = r2 - r1;
		final int dy = -dx - dz;
		return (Math.abs(dx) + Math.abs(dy) + Math.abs(dz))/2;
	}
	
	public int peucd2qr(int q1, int r1, int q2, int r2)
	{
		final int x1 = q1 - (r1 - (r1&1)) / 2;
		final int x2 = q2 - (r2 - (r2&1)) / 2;
		final int dx = x2 - x1;
		final int dz = r2 - r1;
		final int dy = -dx - dz;
		
		final int[] xyz = new int[]{dx,dy,dz};
		periodic_xyz(xyz);
		//a^2+b^2+ab == (a^2 + b^2 + (-a-b)^2)/2
		return (xyz[0]*xyz[0] + xyz[2]*xyz[2] + xyz[0]*xyz[2]);  
	}
	
	public int peucd2xyz(int[] xyzA, int[] xyzB )
	{
		final int[] xyz = new int[]{xyzB[0]-xyzA[0], xyzB[1]-xyzA[1], xyzB[2]-xyzA[2] };
		periodic_xyz(xyz);
		//a^2+b^2+ab == (a^2 + b^2 + (-a-b)^2)/2
		return (xyz[0]*xyz[0] + xyz[2]*xyz[2] + xyz[0]*xyz[2]);  
	}
	
	public int to_qr(final int[] xyz)
	{
		final int dx = xyz[0];
		final int dz = xyz[2];
		final int x = dx + cx;
		final int z = dz + cz;
		// cube -> odd-r
		final int q = x + ( z-(z&1) ) / 2;
		final int r = z;
		return q + r * W;
	}
	
	public int shift_qr(int qr, int[] shift)
	{
		final int q = qr % W;
		final int r = qr / W;
		int x = q - (r - (r&1)) / 2 - cx + shift[0];
		int z = r - cz + shift[2];
		int y = -x - z;
		 
		final int[] xyz = new int[]{x,y,z};
		
		periodic_xyz(xyz);
		
		x = xyz[0] + cx;
		z = xyz[2] + cz;
		
		final int nq = x + ( z-(z&1) ) / 2;
		final int nr = z;
		return nq + nr * W;	
		
		/*
		boolean at_home = true;
		//apply periodic boundary conditions
		do
		{
			at_home = true;
			if(x >= R)
			{
				at_home = false;
				x -= (2*R); 
				y -= (-R);
				z -= (-R);
			}
			else if (x < -R)
			{
				at_home = false;
				x += (2*R);
				y += (-R);
				z += (-R);
			}
		
			if(y > R)
			{
				at_home = false;
				x -= (-R);
				y -= (2*R);
				z -= (-R);
			}
			else if (y <= -R)
			{
				at_home = false;
				x += (-R);
				y += (2*R);
				z += (-R);
			}
		
			if(z >= R)
			{
				x -= (-R);
				y -= (-R);
				z -= (2*R);
			}
			else if (z < -R)
			{
				x += (-R);
				y += (-R);
				z += (2*R);
			}	
		} while (!at_home);
		
		x += cx;
		z += cz;
		
		final int nq = x + ( z-(z&1) ) / 2;
		final int nr = z;
		return nq + nr * W;
		*/	
	}
	
	public int periodic_qr_safe(int q, int r)
	{
		int x = q - (r - (r&1)) / 2 - cx;
		int z = r - cz;
		int y = -x - z;
		
		final int[] xyz = new int[]{x,y,z};
		
		periodic_xyz(xyz);
		
		x = xyz[0] + cx;
		z = xyz[2] + cz;
		
		final int nq = x + ( z-(z&1) ) / 2;
		final int nr = z;
		return nq + nr * W;	
	}
	
	public int periodic_qr(int qr)
	{
		final int q = qr % W;
		final int r = qr / W;
		
		int x = q - (r - (r&1)) / 2 - cx;
		int z = r - cz;
		int y = -x - z;
		
		final int[] xyz = new int[]{x,y,z};
		
		periodic_xyz(xyz);
		
		x = xyz[0] + cx;
		z = xyz[2] + cz;
		
		final int nq = x + ( z-(z&1) ) / 2;
		final int nr = z;
		return nq + nr * W;	
		/*
		boolean at_home = true;
		//apply periodic boundary conditions
		do
		{
			at_home = true;
			if(x >= R)
			{
				at_home = false;
				x -= (2*R); 
				y -= (-R);
				z -= (-R);
			}
			else if (x < -R)
			{
				at_home = false;
				x += (2*R);
				y += (-R);
				z += (-R);
			}
		
			if(y > R)
			{
				at_home = false;
				x -= (-R);
				y -= (2*R);
				z -= (-R);
			}
			else if (y <= -R)
			{
				at_home = false;
				x += (-R);
				y += (2*R);
				z += (-R);
			}
		
			if(z >= R)
			{
				x -= (-R);
				y -= (-R);
				z -= (2*R);
			}
			else if (z < -R)
			{
				x += (-R);
				y += (-R);
				z += (2*R);
			}	
		} while (!at_home);
		
		x += cx;
		z += cz;
		
		final int nq = x + ( z-(z&1) ) / 2;
		final int nr = z;
		return nq + nr * W;
		*/ 
			
	}
	
	public void periodic_xyz(final int[] xyz)
	{
		//for higher dimensions we could iterate and distinguish even and odd indices but there are atm only 3
		boolean at_home = true;
		do
		{
			at_home = true;
			if(xyz[0] >= R)
			{
				at_home = false;
				xyz[0] -= (2*R); 
				xyz[1] -= (-R);
				xyz[2] -= (-R);
			}
			else if (xyz[0] < -R)
			{
				at_home = false;
				xyz[0] += (2*R);
				xyz[1] += (-R);
				xyz[2] += (-R);
			}
			
			if(xyz[1] > R)
			{
				at_home = false;
				xyz[0] -= (-R);
				xyz[1] -= (2*R);
				xyz[2] -= (-R);
			}
			else if (xyz[1] <= -R)
			{
				at_home = false;
				xyz[0] += (-R);
				xyz[1] += (2*R);
				xyz[2] += (-R);
			}
			
			if(xyz[2] >= R)
			{
				xyz[0] -= (-R);
				xyz[1] -= (-R);
				xyz[2] -= (2*R);
			}
			else if (xyz[2] < -R)
			{
				xyz[0] += (-R);
				xyz[1] += (-R);
				xyz[2] += (2*R);
			}	
			
		}
		while(!at_home);
		return;
	}
	
	public boolean is_inside(int[] xyz)
	{
		return  (xyz[0] >= -R && xyz[0] <  R &&
				 xyz[1] >  -R && xyz[1] <= R &&
				 xyz[2] >= -R && xyz[2] <  R );
	}
	
	public boolean on_canvas(int[] xyz)
	{
		final int dx = xyz[0];
		final int dz = xyz[2];
		final int x = dx + cx;
		final int z = dz + cz;
		// cube -> odd-r
		final int q = x + ( z-(z&1) ) / 2;
		final int r = z;
		return ( q >= 0) && (r >= 0) && (q < W) && (r < W);
	}
	
	public int[] to_qr_pair(int[] xyz)
	{
		final int dx = xyz[0];
		final int dz = xyz[2];
		final int x = dx + cx;
		final int z = dz + cz;
		// cube -> odd-r
		final int q = x + ( z-(z&1) ) / 2;
		final int r = z;
		return new int[]{q,r};
	}
	
	public boolean is_inside(int qr)
	{
		return hex_mask[qr];
	}
	
	private boolean inside(int qr)
	{
		final int q = qr % W;
		final int r = qr / W;
		final int x = q - (r - (r&1)) / 2 - cx;
		final int z = r - cz;
		final int y = -x - z;
		return  (x >= -R && x <  R &&
				 y >  -R && y <= R &&
				 z >= -R && z <  R );
	}
}
