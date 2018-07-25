package krami;
//these are needed for the addMark methods of Bond and Ring
import ij.*; //needed for ij.Macro
import ij.gui.*;
import ij.gui.Roi;
import ij.gui.Line;
import ij.gui.Overlay;
import java.awt.*;
import java.awt.Color;
import java.awt.geom.*;
//import java.awt.geom.Line2D.*;
//import java.awt.geom.Line2D.Float;
//needed for managing Vertices,Edges and Faces
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;

public class SurfaceMesh
{
	public static boolean allow_quadrangles_default = false; //more useful for tilted samples
	public static boolean flat_triangle_check_default = false;  
	public static float max_local_dsqr_variation_default = 9.0f;
	public static boolean debug_mode_default = false;
	public static boolean hex_pixels_default = false;
	public static boolean periodic_default = false;
	public static boolean inverted_default = false;
	public static boolean gray_edges_default = true;
	public static int max_neighbors_default = 3;
	
	public boolean allow_quadrangles = allow_quadrangles_default;
	public boolean flat_triangle_check = flat_triangle_check_default;  
	public float max_local_dsqr_variation = max_local_dsqr_variation_default;
	public boolean debug_mode = debug_mode_default;
	public boolean hex_pixels = hex_pixels_default;
	public boolean periodic = periodic_default;
	public boolean inverted = inverted_default;
	public boolean gray_edges = gray_edges_default;
	public int max_neighbors = max_neighbors_default;		
	
	
	public static Color color_atom = Color.BLUE;
	public static Color color_atom_sig = Color.GREEN;
	public static Color color_atom_err = Color.RED.brighter();
	public static Color color_bond = Color.RED;
	public static Color color_spring = Color.GREEN.darker().darker();
	public static Color color_dual_spring = Color.YELLOW.darker();
	
	public static Color color_3ring = Color.RED;
	public static Color color_4ring = Color.PINK;
	public static Color color_5ring = Color.CYAN.darker();
	public static Color color_6ring = Color.YELLOW.darker();
	public static Color color_7ring = Color.MAGENTA;
	public static Color color_8ring = Color.PINK.darker();
	public static Color color_9ring = Color.RED.darker();
	public static Color color_badring = Color.GREEN;
	
	public static Color color_edges = Color.CYAN;
	
	public static SurfaceMesh master_mesh = null;
	private ArrayList<Atom> atoms = new ArrayList<Atom>(1000);
	private ArrayList<Bond> bonds = new ArrayList<Bond>(1500);
	private ArrayList<Ring> rings = new ArrayList<Ring>(300);
	private ArrayList<Spring> springs = new ArrayList<Spring>(1000);
	public Signature signature = null;
	
	public Atom[] fatoms = null;
	public Bond[] fbonds = null;
	public Ring[] frings = null;
	public Spring[] fsprings = null;
	
	public HexPixels imgHP = null;
	
	private boolean needs_update = true; //no atoms but we still need to update the fatoms ...
	
	
	public int interior_bonds = -1;
	public int num_unique_pairs = -1;
	private int impWidth = 512;
	private int impHeight = 512;
	private int impSlice = 1;
	
	public SurfaceMesh(int _impWidth, int _impHeight, int _impSlice)
	{
		this.impWidth = _impWidth;
		this.impHeight = _impHeight;
		this.impSlice = _impSlice;
	}
	
	public void clear()
	{
		incorporate(false);
		
		atoms.clear();
		bonds.clear();
		rings.clear();
		springs.clear();
		
		fatoms = null;
		fbonds = null;
		frings = null;
		fsprings = null;
		
		signature = null;
		needs_update = true;
	}
	
	
	public void clean_up_master_mesh()
	{
		if(master_mesh == null)
		{	return;} //nothing to do
		
		ArrayList<Atom> atoms2 = new ArrayList<Atom>(master_mesh.atoms.size());
		for(int i = 0; i < master_mesh.atoms.size(); ++i)
		{
			final Atom atomi = master_mesh.atoms.get(i);
			if(!atomi.detached)
			{	
				atomi.id = atoms2.size();
				atoms2.add(atomi);	
			}
		}
		master_mesh.atoms = atoms2; //purge all detached atoms
		master_mesh.fatoms = atoms2.toArray( new Atom[0] );
		
		ArrayList<Bond> bonds2 = new ArrayList<Bond>(master_mesh.bonds.size());
		for(int i = 0; i < master_mesh.bonds.size(); ++i)
		{
			final Bond bondi = master_mesh.bonds.get(i);
			if(!bondi.detached)
			{	
				bondi.id = bonds2.size();
				bonds2.add(bondi);	
			}
		}
		master_mesh.bonds = bonds2; //purge all detached bonds
		master_mesh.fbonds = bonds2.toArray( new Bond[0] );
		
		ArrayList<Ring> rings2 = new ArrayList<Ring>(master_mesh.rings.size());
		for(int i = 0; i < master_mesh.rings.size(); ++i)
		{
			final Ring ringi = master_mesh.rings.get(i);
			if(!ringi.detached)
			{	
				ringi.id = rings2.size();
				rings2.add(ringi);	
			}
		}
		master_mesh.rings = rings2; //purge all detached rings
		master_mesh.frings = rings2.toArray( new Ring[0] );
		
		ArrayList<Spring> springs2 = new ArrayList<Spring>(master_mesh.springs.size());
		for(int i = 0; i < master_mesh.springs.size(); ++i)
		{
			final Spring springi = master_mesh.springs.get(i);
			if(!springi.detached)
			{	
				springi.id = springs2.size();
				springs2.add(springi);	
			}
		}
		master_mesh.springs = springs2; //purge all detached springs
		master_mesh.fsprings = springs2.toArray( new Spring[0] );
		
		
		
		needs_update = false;
		
	}
	
	
	//give a positive rating for reasonably unique neighborhoods
	public static int rate_signature(Atom a0)
	{
		final int rs = a0.rings.size();
		if(rs != 3)
		{	return -2;}
		for(int i = 0; i < rs; ++ i)
		{
			if (!a0.rings.get(i).is_interior)
			{	return -1;}
		}
		int r5 = 0;
		int r6 = 0;
		int r7 = 0;
		for(int i = 0; i < rs; ++ i)
		{
			switch ( a0.rings.get(i).edges.size() )
			{
				case 5:	++r5;	break;
				case 6:	++r6;	break;
				case 7:	++r7;	break;
				
			}
		}
		if( (r5 != 1) || (r6 != 1) || (r7 != 1) )
		{
			return 0;	
		}
		return 2;
	}
	
	public static int eval_sig_strings(String sig_str1, String sig_str2)
	{
		String[] ring_defs1 = sig_str1.split("\\n");
		String[] ring_defs2 = sig_str2.split("\\n");
		int score = 0;
		int max_lvl = Math.min(ring_defs1.length,ring_defs2.length);
		for( int lvl = 0; lvl < max_lvl; ++lvl)
		{
			if(ring_defs1[lvl].equals(ring_defs2[lvl]))
			{
				//quite fair approximation for number of matching rings
				score += (1+ring_defs1[lvl].length())/2;
			}
			else
			{	return score;}
		}
		return score;
	}
	
	public int[] ixyz(float[] pos)
	{
		int[] cube = new int[] 
		{
			Math.round(pos[0]),
			0,
			Math.round(pos[2])
		};
		cube[1] = -cube[0]-cube[2];
		return cube;	
	}
	
	
	private Line2D.Float[] get_lines(float[] pos, float[] dest, boolean shorten)
	{
		if(!hex_pixels)
		{
			//apparent distance seems to fit without cutting the periodic frame
			final float dx = dest[0]-pos[0];
			final float dy = dest[1]-pos[1];
			if( (Math.abs(dx) <= impWidth/2) && ((Math.abs(dy) <= impHeight/2)) )
			{
				//for now just ignore periodic non hex images
				return new Line2D.Float[]{new Line2D.Float( pos[0],pos[1], dest[0],dest[1] )}; 
			}
			else
			{
				final float f = shorten?0.5f:1.0f;
				float[] vec = vecTo(pos, dest); //shortest of all possible vectors
				final float pmx = pos[0] + f*vec[0];
				final float pmy = pos[1] + f*vec[1];
				final float dmx = dest[0] - f*vec[0];
				final float dmy = dest[1] - f*vec[1];
				return new Line2D.Float[]
					{	new Line2D.Float( pos[0],pos[1], pmx,pmy ), 
						new Line2D.Float( dest[0],dest[1], dmx,dmy )
					};
			}
		}
		
		
		int[] ip1 = ixyz(pos);
		int[] ip2 = ixyz(dest);
		
		int qr1 = imgHP.to_qr(ip1);
		int q1 = qr1%impWidth;
		int r1 = qr1/impWidth;
		
		int qr2 = imgHP.to_qr(ip2);
		int q2 = qr2%impWidth;
		int r2 = qr2/impWidth;
	
		int[] dp = new int[]{ip2[0]-ip1[0],ip2[1]-ip1[1],ip2[2]-ip1[2]};
		
		if(imgHP.is_inside(dp))
		{
			return new Line2D.Float[]{ new Line2D.Float( q1,r1, q2,r2 )};
		}
		else //resort to only drawing halfs of the line, that are fully on the canvas
		{	
			final float f = shorten?0.5f:1.0f;
			
			int mq1 = q1;
			int mr1 = r1;
			int mq2 = q2;
			int mr2 = r2; 
		
			float[] vec = vecTo(pos,dest); //shortest of all possible vectors
			vec[0] = f*vec[0]+pos[0];
			vec[1] = f*vec[1]+pos[1];
			vec[2] = f*vec[2]+pos[2];	
			int[] pcube = ixyz(vec); //may or may not be on canvas
			
			vec = vecTo(dest,pos); //shortest of all possible vectors
			vec[0] = f*vec[0]+dest[0];
			vec[1] = f*vec[1]+dest[1];
			vec[2] = f*vec[2]+dest[2];	
			int[] dcube = ixyz(vec); //may or may not be on canvas
			
			int[] qrm2 = imgHP.to_qr_pair(pcube);
			mq2 = qrm2[0];
			mr2 = qrm2[1];
			int[] qrm1 = imgHP.to_qr_pair(dcube);
			mq1 = qrm1[0];
			mr1 = qrm1[1];
			Line2D.Float[] lines = new Line2D.Float[2]; 
			lines[0] = new Line2D.Float( q1,r1, mq2,mr2 );
			lines[1] = new Line2D.Float( q2,r2, mq1,mr1 );
			return lines;
			
		}
			
	}
	public Atom getAtom(int id)
	{
		if(id < atoms.size())
		{	return atoms.get(id);}
		else
		{	return null;}
	}
	
	public Atom[] findAtoms(float pos[], float r)
	{
		ArrayList<Atom> candidates = new ArrayList<Atom>(4);
		final float max_dist_sqr = r*r;
		float shortest_distsqr = Float.NaN;
		float longest_distsqr = Float.NaN;
		for(int i = 0; i < atoms.size(); ++i)
		{
			final Atom atomi = atoms.get(i);
			final float dist_sqr = (float)distsqr(pos, atomi.pos);
			if( dist_sqr <= max_dist_sqr)
			{
				if(! (dist_sqr < longest_distsqr) )
				{
					candidates.add(atomi);
					longest_distsqr = dist_sqr;
					if(Float.isNaN(shortest_distsqr) ) 
					{	shortest_distsqr = dist_sqr;}
				}
				else if ( !( dist_sqr > shortest_distsqr) )
				{
					candidates.add(0,atomi);
					shortest_distsqr = dist_sqr;
				}
				else 
				{
					for(int k = 1; k<candidates.size()-1; ++k)
					{
						if(distsqr(pos, candidates.get(k).pos) > dist_sqr)
						{
							candidates.add(k,atomi);
							break;
						}
					}
				}	
			}
		}
		return candidates.toArray(new Atom[0]);
	}
	
	
	public Atom[] get_observers( Atom ma)
	{
		//There ought to be only one, but better to be future proof
		ArrayList<Atom> candidates = new ArrayList<Atom>(4);
		for(int i = 0; i < atoms.size(); ++i)
		{
			final Atom atomi = atoms.get(i);
			if(atomi.master == ma)
			{	candidates.add(atomi);}
		}
		if(debug_mode && candidates.size()>1)
		{
			System.out.println("OUCH, there are multiple obersvers in a single mesh");
		}
		return candidates.toArray(new Atom[0]);
	}
	
	public Atom[] boxed_atoms( float x1, float y1, float x2, float y2)
	{
		//There ought to be only one, but better to be future proof
		ArrayList<Atom> candidates = new ArrayList<Atom>(12);
		for(int i = 0; i < atoms.size(); ++i)
		{
			final Atom atomi = atoms.get(i);
			if(atomi.master != null)
			{	
				float[] pos = atomi.pos;
				if(x1<=pos[0] && pos[0]<=x2 && y1<=pos[1] && pos[1]<=y2)
				{	candidates.add(atomi);}
			}
		}
		return candidates.toArray(new Atom[0]);
	}
	
	
	
	
	public class Atom
	{
		public boolean detached = false;
		public float[] pos = null;
		public int id = -1;
		public boolean is_interior = false;
		public Atom master = null;
		public int observers = 0;
		
		public int ind = 0;
		public float intensity = 1.0f;
		public ArrayList<Atom> neighbors = new ArrayList<Atom>(8);
		public ArrayList<Atom> candidates = new ArrayList<Atom>(34);
		public ArrayList<Bond> edges = new ArrayList<Bond>(8);
		public ArrayList<Ring> rings = new ArrayList<Ring>(8);
		public Atom[] unique_pairs = null;
		public float[] unique_dist_scales = null;
		int top_candidate = -1;
		public float shortest_distsqr = Float.NaN;
		float short_distsqr = Float.NaN;
		float long_distsqr = Float.NaN;
		
		private Atom()
		{
			id = atoms.size();
			atoms.add(this);
		}
		
		
		Bond get_other_bond(Ring r, Bond b)
		{
			if( (!rings.contains(r)) || (!edges.contains(b)) || (!r.has_edge(b)) )
			{	return null;}
			for(int i = 0; i < edges.size(); ++i)
			{
				final Bond bi = edges.get(i);
				if( (bi != b) && (r.has_edge(bi)) && (r.has_vertex(bi.a1)) && (r.has_vertex(bi.a2)) )
				{	return bi;}	
			}
			return null;
		}
		
		
		int addCandidate( Atom other)
		{
			if(other == this)
			{
				System.out.println("Atom cannot be his own candidate");
				return -1;
			}
			
			final float other_distance = distsqrTo(other);
			if(candidates.size()<16 || !(other_distance > long_distsqr))
			{
				int r = 0;
				while( (r < candidates.size()) && (other_distance > distsqrTo(candidates.get(r))) )
				{++r;}
				if(r < 32)
				{
					int or = 0;
					while( (or < other.candidates.size() ) && ( other_distance > other.distsqrTo(other.candidates.get(or))) )
					{++or;}
					if( or < 32)
					{
						candidates.add(r,other);
						other.candidates.add(or,this);
						
						if(candidates.size()>32)
						{
							Atom ca = candidates.remove(candidates.size()-1);
							ca.candidates.remove(this);
							ca.update_dsqrlimits();
								
						}
						
						if(other.candidates.size()>32)
						{
							Atom oca = other.candidates.remove(other.candidates.size()-1);
							oca.candidates.remove(other);
							oca.update_dsqrlimits();	
						}
						update_dsqrlimits();
						other.update_dsqrlimits();
					}
					
				}
			
			 
			}
			/*
			if(candidates.isEmpty())
			{
				System.out.println("Atom: " + id + " near x,y: " + (int)pos[0] + "," + (int)pos[1] + 
				" has no candidates AFTER adding a candidate");
				System.out.println("Candidate: " + other.id + " near x,y: " + (int)other.pos[0] +
				 "," + (int)other.pos[1] + " dist: " + other_distance);
				System.out.println("closest_distance: " + shortest_distance);
				System.out.println("furthest_distance: " + furthest_distance);
				throw new RuntimeException(Macro.MACRO_CANCELED);
			}
			*/
			top_candidate = (candidates.isEmpty()?-1:0);
			
			return candidates.size();
			
		}
		
		void update_dsqrlimits()
		{
			if(!candidates.isEmpty())
			{
				top_candidate = 0;
				short_distsqr = distsqrTo(candidates.get(0));
				long_distsqr = distsqrTo(candidates.get(candidates.size()-1));
			}
			else
			{
				top_candidate = -1;
				short_distsqr = Float.NaN;
				long_distsqr = Float.NaN;
			}	
		}
		
		void addMark(Overlay overlay)
		{
			if(!detached)
			{
				float d = (shortest_distsqr>1) ? (float)Math.sqrt(0.1*shortest_distsqr): -1.0f;
				OvalRoi circle = null;
				float q = pos[0];
				float r = pos[1];
				if(hex_pixels)
				{
					int qr = imgHP.to_qr(ixyz(pos));
					q = qr%impWidth;
					r = qr/impWidth;
				}
				
				
				if(d > 1.0f)
				{	
					circle = new OvalRoi(q-d/2,r-d/2,d,d);
					circle.setStrokeColor( color_atom );	
				}
				else
				{
					d = 6;
					circle = new OvalRoi(q-d/2,r-d/2,d,d);
					circle.setStrokeColor( color_atom_err );	
				}
				
				
				
				if(gray_edges && (!is_interior))
				{	circle.setStrokeColor( Color.LIGHT_GRAY );}	
				if(impSlice != 0)
				{	circle.setPosition(impSlice);}
				if( (signature!=null) && (this == signature.atom) )
				{ circle.setStrokeColor( color_atom_sig );}
				
				if(detached)
				{	circle.setStrokeColor(color_edges);}
				
				
				overlay.add( circle ,"Atom:" + id);
			}
		}
		
		boolean check_interior()
		{
			is_interior = (neighbors.size() > 2);
			for(int i = 0; i < neighbors.size() && is_interior; ++i)
			{	is_interior &= (neighbors.get(i).neighbors.size() > 2);}
			return is_interior;
		}
		
		
		int addNeighbor( Atom other )
		{
			if( ! ( this == other || neighbors.contains(other) || other.neighbors.contains(this) ) )
			{ 
				if( ( (!inverted) && (accept_neighbor(other) && other.accept_neighbor(this) && (neighbors.size() < max_neighbors) ) ) || //rules for atoms
					(   inverted  && (accept_neighbor(other) || (other.accept_neighbor(this)))	) ) //rules for dual atoms aka rings
				{
					new Bond(this, other);
				}
				candidates.remove(other);
				update_dsqrlimits();
				other.candidates.remove(this);
				other.update_dsqrlimits();
						
			}
			return neighbors.size();
		}
		
		float distsqrTo(Atom other)
		{
			return (float)distsqr(pos, other.pos);
		}
		
		
		float distanceTo (Atom other)
		{
			return (float) Math.sqrt(distsqr(pos, other.pos));
		}
		/*
		float distanceTo2D (Atom other)
		{
			return (float) Math.sqrt(dist2Dsqr(pos, other.pos));
		}
		*/
		
		public float[] vectorTo (Atom other)
		{
			/*
			float dx =  other.pos[0] - pos[0];
			float dy =  other.pos[1] - pos[1];
			float dz =  other.pos[2] - pos[2];
			float[] vec = {dx,dy,dz};
			*/
			return vecTo(pos, other.pos);
		}
		
		
		boolean accept_neighbor(Atom other)
		{ 
			if(inverted) return true; //Atoms are actually negative rings, there are no topological restrictions
			if(max_local_dsqr_variation > 0.0)
			{
				if( (short_distsqr > max_local_dsqr_variation * shortest_distsqr  ) ) //magic number
				{	return false;}
			}
			//quick test as more than 4 neighbors would anyways always cause a bond collision
			if(neighbors.size() >= max_neighbors)//no more than 4 bonds
			{	return false;	}
			
			int n2 = other.neighbors.size(); //other is quite not yet a neighbor
			//look for triangles
			for(int n0=0; (n0 < neighbors.size()); ++n0)
			{
				if(neighbors.get(n0).neighbors.contains(other))
				{ return false;}
				n2+=neighbors.get(n0).neighbors.size()-1;
			}
			
			{	//look for quadrangles
				int quads = 0;
				Atom[] secondNeighbors = new Atom[n2];
				int n2w=0;
				for(int n0 = -1; n0 < neighbors.size(); ++n0)
				{
					Atom first = (n0==-1) ? other : neighbors.get(n0);
					for(int n1 = 0; n1 < first.neighbors.size(); ++n1)
					{
						Atom second = first.neighbors.get(n1);
						if(second != this )
						{
							for(int i=0; i < n2w; ++i)
							{
								if(second == secondNeighbors[i])
								{	++quads;}	
							}
							secondNeighbors[n2w++] = second;
						}
					}
				}
				if(quads > 2) //two quadrangles mean that a hexagon was divided
				{	return false;}
				else if(quads > 0 && !allow_quadrangles)
				{	return false;}
			}
			if(flat_triangle_check && neighbors.size()==2) //Check if the third and first two neighbors would form a 2D triangle around us
			{
				// Compute vectors        
				float[] v0 = neighbors.get(0).vectorTo(other);
				float[] v1 = neighbors.get(0).vectorTo(neighbors.get(1));
				float[] v2 = neighbors.get(0).vectorTo(this);
				
				// Compute dot products				
				final float dot00 = v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2];
				final float dot01 = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
				final float dot02 = v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2];
				final float dot11 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
				final float dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
				
				// Compute barycentric coordinates
				final float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
				final float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
				final float v = (dot00 * dot12 - dot01 * dot02) * invDenom;
				
				// Check if we are inside triangle, or the tip of the pyramid hover above the ground
				return (u > 0.0f) && (v > 0.0f) && (u + v < 1.0f);	
			}
			
			
				
			return true;
		}
		
		double get_roughness()
		{
			if(neighbors.size()==3)
			{
				float[] u = vectorTo(neighbors.get(0));
				float[] v = vectorTo(neighbors.get(1));
				float[] w = vectorTo(neighbors.get(2));
				
				double[] cr = new double[3];
				cr[0] = u[1]*v[2] - u[2]*v[1];
				cr[1] = u[2]*v[0] - u[0]*v[2];
				cr[2] = u[0]*v[1] - u[1]*v[0];  
				
				return Math.abs( w[0]*cr[0] + w[1]*cr[1] + w[2]*cr[2] );
			}
			else
			{	return 0.0;}
			
		}
		
		double get_heronShape()
		{
			if(neighbors.size() != 3)
			{	return -1.0;}
			final double a = neighbors.get(0).distanceTo(neighbors.get(1));
			final double b = neighbors.get(0).distanceTo(neighbors.get(2));
			final double c = neighbors.get(1).distanceTo(neighbors.get(2));
			
			final double is = 2/(a+b+c);
			return (1.0-a*is)*(1.0-b*is)*(1.0-c*is);	
		}
		
		int put_unique_distances( float[] global_distances, int nd)
		{
			for(int n0 = 0; n0 < unique_pairs.length; ++n0)
			{
				Atom other = unique_pairs[n0];
				global_distances[nd++] = distanceTo(other)*unique_dist_scales[n0];
			}	
			return nd;
		}
		
		int init_unique_pairs()
		{
			int nd=0;
			//first run just count how many there are
			for(int n0 = 0; n0 < neighbors.size(); ++n0)
			{
				Atom first = neighbors.get(n0);
				if(ind < first.ind)
				{
					nd++;
				}
				for(int n1 = 0; n1 < first.neighbors.size(); ++n1)
				{
					Atom second = first.neighbors.get(n1);
					if(ind < second.ind )
					{
						nd++;
					}
				}
			}	
			
			unique_pairs = new Atom[nd];
			unique_dist_scales = new float[nd];
			//second run fill in the arrays
			nd=0;
			for(int n0 = 0; n0 < neighbors.size(); ++n0)
			{
				Atom first = neighbors.get(n0);
				if(ind < first.ind)
				{
					unique_dist_scales[nd] = 1.0f;
					unique_pairs[nd++] = first;
				}
				for(int n1 = 0; n1 < first.neighbors.size(); ++n1)
				{
					Atom second = first.neighbors.get(n1);
					if(ind < second.ind )
					{
						unique_dist_scales[nd] = (float)(1.0/Math.sqrt(3.0));
						unique_pairs[nd++] = second;
					}
				}
			}	
			return nd;
		}
		
		void sort_neighbors()
		{
			//simple sort suffice for for N <= 8
			//sort the edges along with the neighbors as they should
			//highly correlated with the neighbors
			for(int i = 0; i<neighbors.size(); ++i)
			{	
				float[] veci = vectorTo(neighbors.get(i));
				double thetai = Math.atan2(-veci[1],veci[0]);
				for(int j = i+1; j<neighbors.size(); ++j)
				{
					float[] vecj = vectorTo(neighbors.get(j));
					final double thetaj = Math.atan2(-vecj[1],vecj[0]);
					if(thetaj < thetai)
					{
						Collections.swap(neighbors,i,j);
						Collections.swap(edges,i,j);
						thetai = thetaj;
					}
					
				}
			}
			
			for(int i = 0; i<edges.size(); ++i)
			{	
				float[] veci = vecTo(pos,edges.get(i).pos);
				double thetai = Math.atan2(-veci[1],veci[0]);
				for(int j = i+1; j<edges.size(); ++j)
				{
					float[] vecj = vecTo(pos,edges.get(j).pos);
					final double thetaj = Math.atan2(-vecj[1],vecj[0]);
					if(thetaj < thetai)
					{
						Collections.swap(edges,i,j);
						thetai = thetaj;
					}
					
				}
			}
			
			
		}
		
		public void detach()
		{
			if(!detached)
			{
				final int es = edges.size();
				for(int i = es-1; i >= 0; --i)
				{	edges.get(i).detach();}
				id = -1;
				intensity = 0.0f;
				detached = true;
				needs_update = true;
			}
		}
		
		
		private Atom ( float[] xyz)
		{
			pos = xyz;
			
			if(hex_pixels)
			{	ind = imgHP.to_qr(ixyz(pos)); }
			else
			{	ind = Math.round(pos[0]) + impWidth * Math.round(pos[1]);}
			needs_update = true;
			id = atoms.size(); 
			atoms.add(this);	
		}
		
		public void debug_info()
		{
			System.out.println("DEBUG atom: " + id + " at: (" + (int)pos[0] +  "," + (int)pos[1] + ")" +
			( (impSlice!=0)?(" in mesh: " + impSlice): " in master_mesh" ) );
			System.out.print("touching rings: [ ");
			final int ns = rings.size();
			for(int i = 0; i < ns; ++i)
			{	System.out.print( (i==0?"": "-") + rings.get(i).id);}
			System.out.println("]");
			System.out.print("neighbors: [ ");
			final int vs = neighbors.size();
			for(int i = 0; i < vs; ++i)
			{	System.out.print( (i==0?"": ", ") + i + "(" + (int)neighbors.get(i).pos[0]+ "," + (int)neighbors.get(i).pos[1] + ")"  );}
			System.out.println("]");
			System.out.print("bonds: [ ");
			final int es = edges.size();
			for(int i = 0; i < es; ++i)
			{	
				final Atom a1 = neighbors.get(i);
				final Atom a2 = edges.get(i).get_other_atom(this);
				System.out.print( ((i==0)?"":", ") + "to " + i + ((a1==a2)?"(ok)":"(broken)")  );	
			}
			System.out.println("]");
		}
		
	}
	
	public class Bond
	{
		public boolean detached = false;
		public float[] pos = null;
		public int id = -1;
		public boolean is_interior = false;
		public Bond master = null;
		public int observers = 0;
		
		public Atom a1 = null;
		public Atom a2 = null;
		public Ring left_ring = null;
		public Ring right_ring = null;
		
		public Spring spring = null;
		private Line2D.Float[] segments = null;
		
		private Bond()
		{
			id = bonds.size();
			bonds.add(this);	
		}
		
		float getLength()
		{
			/*
			final float dx = a2.pos[0] - a1.pos[0];
			final float dy = a2.pos[1] - a1.pos[1];
			final float dz = a2.pos[2] - a1.pos[2];
			return (float) Math.sqrt(dx*dx+dy*dy+dz*dz);
			*/
			return (float)Math.sqrt(distsqr(a1.pos,a2.pos));
		}
		/*
		float getLength2D()
		{
			final float dx = a2.pos[0] - a1.pos[0];
			final float dy = a2.pos[1] - a1.pos[1];
			return (float) Math.sqrt(dx*dx+dy*dy);
		}
		*/
		Atom get_other_atom(Atom a)
		{
			if(a == a1)
			{	return a2;}
			else if(a == a2)
			{	return a1;}
			else
			{	return null;}
			
		}
		
		
		public Ring get_other_ring(Ring r1)
		{
			if(r1 == left_ring)
			{	return right_ring;}
			else if(r1 == right_ring)
			{	return left_ring;}
			else
			{	return null;}
			
		}
		/*
		int getcenter()
		{
			final int cx = (int) (0.5*(a1.pos[0] + a2.pos[0]+1));
			final int cy = (int) (0.5*(a1.pos[1] + a2.pos[1]+1));
			return cx + cy * impWidth;
		}
		*/
		boolean check_interior()
		{
			is_interior = (a1.neighbors.size()>2) && (a2.neighbors.size()>2);
			return is_interior;	
		}
		
		void addMark(Overlay overlay)
		{
			if(!detached)
			{
				
				/*
				if(hex_pixels)
				{
					int[] ip1 = ixyz(a1.pos);
					int[] ip2 = ixyz(a2.pos);
					
					int qr1 = imgHP.to_qr(ip1);
					q1 = qr1%impWidth;
					r1 = qr1/impWidth;
					
					int qr2 = imgHP.to_qr(ip2);
					q2 = qr2%impWidth;
					r2 = qr2/impWidth;
					
					int[] dp = new int[]{ip2[0]-ip1[0],ip2[1]-ip1[1],ip2[2]-ip1[2]};
					if(!imgHP.is_inside(dp)) //resort to only drawing the longer half of the line
					{	
					
						int[] ipm = ixyz(pos);
						int qrm = imgHP.to_qr(ipm);
					
						int[] dp1 = new int[]{ ip1[0]-ipm[0], ip1[1]-ipm[1], ip1[2]-ipm[2]};
						int[] dp2 = new int[]{ ip2[0]-ipm[0], ip2[1]-ipm[1], ip2[2]-ipm[2]};
						if( ( (dp1[0]*dp1[0] + dp1[2]*dp1[2] + dp1[0]*dp1[2]) < 
							(dp2[0]*dp2[0] + dp2[2]*dp2[2] + dp2[0]*dp2[2]) ) )
						{//a1 appears closer to periodic center
							q2 = qrm%impWidth;
							r2 = qrm/impWidth;
						}
						else //just go with a2
						{
							q1 = qrm%impWidth;
							r1 = qrm/impWidth;
						}	
					
					}
				}
				*/
				Line2D.Float[] lines = get_lines(a1.pos,a2.pos, true);
				if(lines != null)
				{
					for(int i = 0; i < lines.length; ++i)
					{
						Line line = new Line ( lines[i].x1, lines[i].y1, lines[i].x2, lines[i].y2 );
						//line.setStrokeColor(is_interior? Color.YELLOW : Color.LIGHT_GRAY);
						line.setStrokeColor(left_ring != null && right_ring != null? color_bond : color_edges);
						if(impSlice != 0)
						{	line.setPosition(impSlice);}
						overlay.add( line,"Bond:" + id);
					}
				}							
			}
		}
		
		boolean collision()
		{
			/*if(inverted) //bonds are actually springs
			{
				return intersecting();
			}*/
			//else //perform a more restrictive exclusive disk check
			//{
				
				final double r2min = distsqr(a1.pos,a2.pos)*0.25d;
				final int bs = bonds.size();
				for(int i = 0; (i < bs); ++i)
				{	
					if(i == id) {continue;}
					final Bond bi = bonds.get(i);
					if(bi.detached) {continue;}
					final boolean shared_atom = (bi.a1 == a1 || bi.a2 == a2 || bi.a1 == a2 || bi.a2 == a1);
					if(inverted && shared_atom )
					{	continue;} //exclude those other "springs" who share an atom with this "spring"
					if(inverted)
					{
						if(intersecting(bi))
						{ 
							double a = getLength();
							double b = bi.getLength();
							/*
							if(Math.abs(a-b)/(a+b) < 0.2) //roughly a square
							{
								bi.detach();
								return true;
							}
							else //remove or cancel the longer of the two
							*/
							
							if(a > b)
							{	return true;}
							else if(a < b)
							{	bi.detach();}
							
						}
						continue;
					}
					
					
					
					
					if (r2min > distsqr(pos,bi.pos) ) 
					{
						return true;		
					}	
				}	
				return false;
			//}
		}
		
		boolean intersecting(Bond bi) 	
		{	
			final Line2D.Float[] msegs = segments;
			
			if(!bi.detached)
			{
				if(inverted && (bi.a1 == a1 || bi.a2 == a2 || bi.a1 == a2 || bi.a2 == a1) )
				{	return false;}
				final Line2D.Float[] osegs = bi.segments;	
				for(int m = 0; m < msegs.length; ++m)
				{
					final Line2D.Float ms = msegs[m];
					for(int o = 0; o < osegs.length; ++o)
					{
						if( ms.intersectsLine(osegs[o]))
						{
							return true;	
						}
					}
				}
			}	
			return false;
		}
		
		
		
		public void detach()
		{
			if(! detached)
			{
				a1.edges.remove(this);
				a2.edges.remove(this);
				a1.neighbors.remove(a2);
				a2.neighbors.remove(a1);
				if(left_ring != null)
				{	left_ring.detach();}
				if(right_ring != null)
				{	right_ring.detach();}
				//a1 = null;
				//a2 = null;
				
				id = -1;
				left_ring = null;
				right_ring = null;
				if(spring != null)
				{
					spring.bond = null;
					spring.detached = true;
				}
				detached = true;
			}	
		}
		
		
		
		private Bond( Atom a, Atom b )
		{
			a1 = a;
			a2 = b;
			pos = vecTo(a1.pos,a2.pos);
			pos[0] = 0.5f*pos[0]+a1.pos[0];
			pos[1] = 0.5f*pos[1]+a1.pos[1];
			pos[2] = 0.5f*pos[2]+a1.pos[2];
			segments = get_lines(a1.pos,a2.pos,false);
			if(hex_pixels)
			{	
				int[] cube = ixyz(pos); 
				imgHP.periodic_xyz(cube);
				pos[0] = (float)cube[0];
				pos[1] = (float)cube[1];
				pos[2] = (float)cube[2];
			}
			if(!collision())
			{
				
				a1.neighbors.add(a2);
				a2.neighbors.add(a1);
				a1.edges.add(this);
				a2.edges.add(this);
				id = bonds.size();
				bonds.add(this);
				
				final float bond_dsqr = (float)distsqr(a1.pos,a2.pos);
				if(! (a1.shortest_distsqr > bond_dsqr) )
				{	a1.shortest_distsqr = bond_dsqr;}
				if(! (a2.shortest_distsqr > bond_dsqr) )
				{	a2.shortest_distsqr = bond_dsqr;}	
			}
			
			
			
		}
		
		public void debug_info()
		{
			System.out.println("DEBUG bond: " + id + ( (impSlice!=0)?(" in mesh: " + impSlice + " master: " + ((master==null)?"null":master.id)): " in master_mesh" ) +
			" at: " + (int)pos[0] + "," + (int)pos[1] +
			 " spring: " + (spring==null?"null":spring.id) );
			System.out.print("touching rings: [ ");
			System.out.print((left_ring!=null)?"left: " + left_ring.id + ", " : "left: null, " );
			System.out.print((right_ring!=null)?"right: " + right_ring.id : "right: null" );
			System.out.println(" ]");
			System.out.print("atoms: [ ");
			System.out.print( (a1!=null)?"a1: " + a1.id +
			  "(" + (int)a1.pos[0]+ "," + (int)a1.pos[1] + ") ":"a1: null ");
			System.out.print( (a1.master!=null)?"master: " + a1.master.id +
			  "(" + (int)a1.master.pos[0]+ "," + (int)a1.master.pos[1] + "), ":"master: null, ");
			System.out.print( (a2!=null)?"a2: " + a2.id +
			  "(" + (int)a2.pos[0]+ "," + (int)a2.pos[1] + ") ":"a2: null ");
			System.out.print( (a2.master!=null)?"master: " + a2.master.id +
			  "(" + (int)a2.master.pos[0]+ "," + (int)a2.master.pos[1] + ") ":"master: null ");
			System.out.println(" ]");
		}
		
	}
	
	public class Ring
	{
		public boolean detached = false;
		public float[] pos = null;
		public int id = -1;
		public boolean is_interior = false;
		public Ring master = null;
		public int observers = 0;
		public float intensity = 1.0f;
		
		public double radius = 0.0;
		public ArrayList<Atom> vertices = new ArrayList<Atom>(10);
		public ArrayList<Bond> edges = new ArrayList<Bond>(10);
		public ArrayList<Ring> neighbors = new ArrayList<Ring>(10);
		public ArrayList<Spring> springs = new ArrayList<Spring>(10);
		public boolean is_closed = false;
		
		private Ring()
		{
			id = rings.size();
			rings.add(this);
		}
	
		public boolean has_edge(Bond bond)
		{	return edges.contains(bond);}
		
		public boolean has_vertex(Atom atom)
		{	return vertices.contains(atom);}
	
		public boolean shared_edge(Ring other)
		{
			if(other==null || edges.isEmpty() || other.edges.isEmpty())//trivial cases
			{	return false;}
			return !Collections.disjoint(edges,other.edges);
		}
		//any common edge should do for the elimination check
		private Bond common_edge(Ring other)
		{
			//if(neighbors.contains(other))
			//{return null;}
			for(int i = 0; i < edges.size(); ++i)
			{
				final Bond bondi = edges.get(i);
				if(bondi.left_ring == other || bondi.right_ring == other)
				{	return bondi;}
			}
			return null;
		}
		
		
	
		public boolean check_interior()
		{
			is_interior = ( edges.size() == neighbors.size());
			/*
			//quite restrictive
			is_interior = ( edges.size() == neighbors.size() );
			//extra restrictive
			for(int n = 0; is_interior && n < neighbors.size(); ++n) 
			{	is_interior = neighbors.get(n).edges.size() == neighbors.get(n).neighbors.size();}
			*/
			return is_interior;
		}
		
		void addMark(Overlay overlay)
		{
			if(!detached)
			{
				double dsqr = 0.0;
				final int es = edges.size();
				for(int i = 0; i < es; ++i)
				{
					dsqr += ( hex_pixels ? distsqr(pos,edges.get(i).pos) : dist2Dsqr(pos,edges.get(i).pos) );
				}
				final float d = (float)(Math.sqrt(1.5*dsqr/es));
				
				float q = pos[0];
				float r = pos[1];
				if(hex_pixels)
				{
					int qr = imgHP.to_qr(ixyz(pos));
					q = qr%impWidth;
					r = qr/impWidth;
				}
				
				
				OvalRoi circle = new OvalRoi(q-d/2,r-d/2,d,d);
				if( is_interior || (!gray_edges) )
				{
					switch(es)
					{
						case 3:
							circle.setStrokeColor( color_3ring );
							break;
						case 4:
							circle.setStrokeColor( color_4ring );
							break;
						case 5: 
							circle.setStrokeColor( color_5ring );
							break;
						case 6:
							circle.setStrokeColor( color_6ring );		
							break;
						case 7: 
							circle.setStrokeColor( color_7ring );
							break;
						case 8:
							circle.setStrokeColor( color_8ring );
							break;
						case 9:
							circle.setStrokeColor( color_9ring );
							break;	
						default:
							circle.setStrokeColor( color_badring );	
					}
				}
				else
				{
					circle.setStrokeColor( Color.LIGHT_GRAY );	
				}
				if(impSlice != 0)
				{	circle.setPosition(impSlice);}
				overlay.add( circle , "Ring:" + id);							
				/*
				for(int i = 0; i < neighbors.size(); ++i)
				{
					Ring other = neighbors.get(i);
					if(id < other.id)
					{	continue;}
					
					Line2D.Float[] lines = get_lines(pos,other.pos, true);
					if(lines != null)
					{
						for(int k = 0; k < lines.length; ++k)
						{
							Line line = new Line ( lines[k].x1, lines[k].y1, lines[k].x2, lines[k].y2 );
							line.setStrokeColor( (is_interior || other.is_interior) ? Color.GREEN.darker().darker() : Color.GRAY );
							if(impSlice != 0)
							{	line.setPosition(impSlice);}
							overlay.add( line );
						}
					}				
				}
				*/
			}
		}
		
		public void detach() //do we need saftey checks here!?
		{
			if(!detached)
			{
				//springs would not yet exist if close() fails
				final int ss = springs.size();
				for(int i = 0; i < ss; ++i)
				{
					springs.get(i).bond.spring = null;
					springs.get(i).bond = null;
					springs.get(i).left_atom = null;
					springs.get(i).right_atom = null;
					Ring other = springs.get(i).get_other_ring(this);
					if(other != null) //should always be true
					{	other.springs.remove(springs.get(i));}
					springs.get(i).detached = true;
				}
				springs.clear();
				final int ns = neighbors.size();
				for(int i = 0; i < ns; ++i)
				{	neighbors.get(i).neighbors.remove(this);}
				neighbors.clear();
				final int es = edges.size();
				for(int i = 0; i < es; ++i)
				{	
					if(edges.get(i).left_ring == this)
					{	edges.get(i).left_ring = null;}
					else if (edges.get(i).right_ring == this)
					{	edges.get(i).right_ring = null;}
				}
				edges.clear();
				final int vs = vertices.size();
				for(int i = 0; i < vs; ++i)
				{	vertices.get(i).rings.remove(this);}	
				vertices.clear();
				detached = true;
				
			}
		}
	
		public void debug_info()
		{
			System.out.println("DEBUG ring: " + id + ( (impSlice!=0)?(" in mesh: " + impSlice): " in master_mesh" ) +
			 " at: " + (int)pos[0] + "," + (int)pos[1] + " closed: " + is_closed);
			/* 
			System.out.print("thetas: [ ");
			final int vs = vertices.size();
			for(int i = 0; i < vs; ++i)
			{	
				float[] veci = vecTo(pos,vertices.get(i).pos);
				double thetai = Math.atan2(-veci[1],veci[0]);	
				System.out.print( (i==0?"": " | ") + (Math.round(100.0 * thetai)/100.0) );
			}
			System.out.println("]");
			*/
			
			final int vs = vertices.size();
			System.out.print("" + vs + " vertices: [ ");
			for(int i = 0; i < vs; ++i)
			{	System.out.print( (i==0?"": "-") + vertices.get(i).id + "(" + (int)vertices.get(i).pos[0]+ "," + (int)vertices.get(i).pos[1] + ")"  );}
			System.out.println("]");
			
			final int es = edges.size();
			System.out.print("" + es + " edges: [ ");
			for(int i = 0; i < es; ++i)
			{	
				final Atom a1 = edges.get(i).a1;
				final Atom a2 = edges.get(i).a2;
				final int idb = edges.get(i).id;
				System.out.print( (i==0?"": "-") + idb + "(" + a1.id + "->" + a2.id + ")"  );
				
			}
			System.out.println("]");
			
			final int ns = neighbors.size();
			System.out.print("" + ns + " neighbors: [ ");
			for(int i = 0; i < ns; ++i)
			{	System.out.print( (i==0?"": "-") + neighbors.get(i).id + "(" + (int)neighbors.get(i).pos[0]+ "," + (int)neighbors.get(i).pos[1] + ")"  );}
			System.out.println("]");
			
			final int ss = springs.size();
			System.out.print("" + ss + " springs: [ ");
			for(int i = 0; i < ss; ++i)
			{	
				final Ring r1 = springs.get(i).r1;
				final Ring r2 = springs.get(i).r2;
				final int idb = springs.get(i).id;
				System.out.print( (i==0?"": "-") + idb + "(" + r1.id + "->" + r2.id + ")"  );
				
			}
			System.out.println("]");
			
			
		}
		
		public boolean add_edge(Bond nb)
		{
			boolean has_tail = false;
			if(is_closed)
			{
				System.out.println("caught adding bond to closed ring in add_edge to ring: " + rings.size());
				return false;
			}
			if( (nb.left_ring != null) && (nb.right_ring != null) )
			{
				System.out.println("caught adding left and right saturrated bond in add_edge to ring: " + rings.size());
				return false;
			}
			if( (nb.left_ring == this) || (nb.right_ring == this) )
			{
				System.out.println("caught adding already marked bond in add_edge to ring: " + rings.size());
				return false;
			}
			if(has_edge(nb))
			{
				System.out.println("caught adding a bond twice (MARKING PASSED!) in add_edge to ring: " + rings.size());
				return false;
			}
			
			int matching_atoms = 0;
			int matching_vertices = 0;
			for(int i = 0; i < edges.size(); ++i)
			{
				final Bond bondi = edges.get(i);
				if(Collections.frequency(edges,bondi) != 1)
				{
					System.out.println("An edge is contained more often than once in ring: " + rings.size());
					return false;
				}
				
				if(bondi == nb)
				{
					System.out.println("caught adding of same bond again (has passed 2 earlier checks!!!)");
					return false;
				}
				if(bondi.a1 == nb.a1 || bondi.a1 == nb.a2)
				{	++matching_atoms;}
				if( bondi.a2 == nb.a1 || bondi.a2 == nb.a2)
				{	++matching_atoms;}
				if( !has_vertex(bondi.a1) || !has_vertex(bondi.a2) )
				{
					System.out.println("an edges ends are no vertices!!!");
					return false;
				}
				if( (Collections.frequency(vertices,bondi.a1) != 1) ||
					(Collections.frequency(vertices,bondi.a2) != 1) )
				{
					System.out.println("An end of an edge appears more than once in vertices");
					return false;
				}	
				
			}
			if(has_vertex(nb.a1))
			{	++matching_vertices;}
			if(has_vertex(nb.a2))
			{	++matching_vertices;}
			if(matching_vertices < matching_atoms)
			{
				//Here we are not just closing but crossing our previous path
				has_tail = true;
				matching_atoms = matching_vertices;
			}
			else if(matching_vertices > matching_atoms)
			{
				System.out.println("Error while adding edge to ring: " + rings.size() +
				 " matching_vertices: " + matching_vertices + " > " + "matching_atoms: " + matching_atoms );
			}
			switch(matching_atoms)
			{
				case 0:
					System.out.println("caught adding of not connected bond");
					return false;
				case 2:
					is_closed = true;
				//falltrough
				case 1:
					if(nb.left_ring == null)
					{	nb.left_ring = this;}
					else if(nb.right_ring == null)
					{	nb.right_ring = this;}
					else
					{
						System.out.println("ERROR, caught adding of double winged bond only AFTER first check");
						return false;
					}
					if(!edges.contains(nb))
					{	
						edges.add(nb);
						Ring other = nb.get_other_ring(this);
						/*
						if(nb.id == 1586)
						{
							System.out.println("adding bond: " + nb.id + " to ring: " + rings.size() + " other ring: " + 
												((	other == null)?"null":other.id) );
						
							if( (other != null) && (other != this) )
							{	
								System.out.println("(other != null) && (other != this) == true");
								if( !neighbors.contains(other) )
								{
									System.out.println("!neighbors.contains(other) == true");
								}
								else
								{
									System.out.println("neighbors size: " + neighbors.size() +
											" do already comprise the other ring: " + other.id);
									System.out.println("************************");	
									debug_info();
									other.debug_info();
									System.out.println("************************");	
								} 
							
							}
						}
						*/
						if( (other != null) && (other != this) )
						{
							if( !neighbors.contains(other) ) //there may be vertices with less than 3 edges
							{	
								neighbors.add(other);
								/*
								if(nb.id == 1586)
								{	System.out.println("confirmed adding neighbor: " + other.id + " total neighbors: " + neighbors.size());}
								*/ 
							}
							else if( (!is_closed) && (nb.a1.edges.size()>2) && (nb.a2.edges.size()>2) )
							{
								System.out.println("Warning, a freshly added edge without any primitive vertex does connect to an already known other ring(1) ");	
							}
							if(!other.neighbors.contains(this))
							{	
								other.neighbors.add(this);
								/*
								if(nb.id == 1586)
								{	System.out.println("confirmed adding reverse neighbor: " + rings.size() + " other total neighbors: " + other.neighbors.size());}
								*/
							}
							else if( (!is_closed) && (nb.a1.edges.size()>2) && (nb.a2.edges.size()>2) )
							{
								System.out.println("Warning, a freshly added edge without any primitive vertex does connect to an already known other ring(2)");	
							}
						
							final double rrdistsqr = distsqr(pos,other.pos);
							final double rbdistsqr = distsqr(pos,nb.pos); 
							if(rrdistsqr/rbdistsqr <= 1.125)
							{
								System.out.println("WARNING, caught adding a bond that is further away than the other ring, aka enclosed");
							}
						
						}	
					}
					else
					{	
						System.out.println("ERROR, caught adding of same bond twice, only AFTER second check");
						return false;			
					}
					
					break;
				default:
					System.out.println("error: a new edge cannot have " + matching_atoms + " matching vertices");
					return false;		
			}
			if(!is_closed)
			{
				if(!has_vertex(nb.a1))
				{
					vertices.add(nb.a1);
					nb.a1.rings.add(this);
					if((nb.a1.edges.size() > 1 ) && (nb.a1.rings.size()>nb.a1.edges.size()))
					{
						System.out.println("Messed up topology of nb.a1 in add_edge to ring: " + rings.size() );
					}	
				}
				else if(!has_vertex(nb.a2))
				{
					vertices.add(nb.a2);
					nb.a2.rings.add(this);
					if((nb.a2.edges.size()>1) && (nb.a2.rings.size()>nb.a2.edges.size()) )
					{
						System.out.println("Messed up topology of nb.a2 in add_edge to ring: " + rings.size() );
					}		
				}
				
			}
			update_pos();
			if(has_tail)
			{	return cut_tail();}
			return true;
		}
		
		boolean cut_tail()
		{
			if(!is_closed)
			{
				System.out.println("Cannot cut tails on open ring: " + rings.size());
				return false;
			}
			boolean eliminate = false;
			do //optional 2nd elimantion pass
			{
				int dropped = 0;
				do
				{
					dropped = 0;
					for(int i=0; i < vertices.size(); ++i)
					{
						int connections = 0;
						Atom atomi = vertices.get(i);
						for(int j = 0; j < atomi.edges.size(); ++j)
						{	
							if( has_edge(atomi.edges.get(j)) )
							{
								++connections;
							}	
						}
						if( eliminate || (connections == 1) )
						{
							vertices.remove(atomi);
							dropped = 1;
							atomi.rings.remove(this);
							
							
							for(int j = 0; j < atomi.edges.size(); ++j)
							{	
								Bond bondj = atomi.edges.get(j);
								if( has_edge(bondj) )
								{
									if(bondj.left_ring == this)
									{	bondj.left_ring = null; }
									if(bondj.right_ring == this)
									{	bondj.right_ring = null;	}	
									edges.remove(bondj);
									if(bondj.left_ring != null)
									{
										if(!shared_edge(bondj.left_ring))
										{
											neighbors.remove(bondj.left_ring);
											bondj.left_ring.neighbors.remove(this);
										}
									}
									if(bondj.right_ring != null)
									{
										if(!shared_edge(bondj.right_ring))
										{
											neighbors.remove(bondj.right_ring);
											bondj.right_ring.neighbors.remove(this);
										}
									}
									
								}
							}
							
							update_pos();
							break; // for vertices was invalidated
						}	
						
					}	
				}
				while(dropped > 0);
				
				if(eliminate) break;
				
				final int ns = neighbors.size();
				for(int i = 0; (i < ns) && (!eliminate); ++i)
				{
					final Ring ringi = neighbors.get(i);
					final float[] rpos = ringi.pos;
					final Bond bondi = common_edge(ringi);
					final float[] bpos = bondi.pos;
					
					final float[] l = vecTo(bondi.a1.pos,bondi.a2.pos);
					final float[] a = vecTo(bpos,pos);
					final float[] b = vecTo(bpos,rpos);
					
					if((a[0]*b[0]+a[1]*b[1]+a[2]*b[2]) > 0.0f)
					{
						eliminate = true;
					}
					
					final double dsqrR = distsqr(pos,rpos);
					final double dsqrB = distsqr(pos,bpos);
					
					if( dsqrR/dsqrB < 1.125) //maybe we could increse limit
					{
						eliminate = true;
					}	
					
				
				
				
				}
				
				///Weiredly enough this still shows up ?
				/*
				if(eliminate)
				{
					System.out.println("suspicous ring id: " + rings.size() + " edges: " + edges.size() +
							 " vertices: " + vertices.size() + " neighbors: " + neighbors.size() + " springs: " + springs.size() +
							  " at x,y: " + (int)pos[0] + "," + (int)pos[1] );
							  debug_info();
					//eliminate = false;
				}
				*/ 
				
				
			}
			while(eliminate);
			return !eliminate;
		}
				
		void update_pos()
		{			
			if(vertices.isEmpty())
			{
				pos = new float[]{0.0f,0.0f,0.0f};
				return;
			}
			Atom atom0 = vertices.get(0);
			float cx = 0.0f;
			float cy = 0.0f;
			float cz = 0.0f;
			final int s = edges.size();
			for(int i=0;i<s;++i)
			{
				Bond bondi = edges.get(i);
				float[] vec = vecTo(atom0.pos,bondi.pos);
				cx += vec[0];
				cy += vec[1];
				cz += vec[2];
			}
			if(s > 1)
			{
				final float is = 1.0f/s;
				cx *= is;
				cy *= is;
				cz *= is;
			}
			pos[0] = atom0.pos[0] + cx;
			pos[1] = atom0.pos[1] + cy;
			pos[2] = atom0.pos[2] + cz;
			if(hex_pixels)
			{
				int[] cube = ixyz(pos);
				imgHP.periodic_xyz(cube);
				pos[0] = (float)cube[0];
				pos[1] = (float)cube[1];
				pos[2] = -pos[0]-pos[1];
			}
			else if(periodic)
			{
				pos[0]=(pos[0] + impWidth)%impWidth;
				pos[1]=(pos[1] + impHeight)%impHeight;
			}
			
		}
		
		boolean close()
		{
			if(is_closed)
			{	return true;}	
			while( !is_closed)
			{
				final int nv = vertices.size();
				int iend = -1;
				for(int j = 0; j < 2; ++j)
				{
					final int i = ( (j==0)? nv-1 : 0);	
					
					Atom atomi = vertices.get(i);
					int shared_edges = 0;
					int free_edges = 0;
					int single_edges = 0;
					int double_edges = 0;
					int closing_edges = 0;
					if( (atomi.edges.size()>1) && (atomi.rings.size() > atomi.edges.size()) )
					{
						System.out.println("Faulty atom near x,y: " + (int)atomi.pos[0] + "," + (int)atomi.pos[1] + 
						" rings: " + atomi.rings.size() + " > edges: " + atomi.edges.size() + " in closing ring: " + rings.size() );
						for(int l=0; l < atomi.rings.size(); ++l)
						{
							System.out.println("SUSPECT: " + l);
							atomi.rings.get(l).debug_info(); 
						}
						throw new RuntimeException(Macro.MACRO_CANCELED);
					}
					if(atomi.edges.size() != atomi.neighbors.size())
					{
						System.out.println("Faulty atom rings: " + atomi.edges.size() + "!=" + atomi.neighbors.size() + " in closing ring: " + rings.size() );
					}
					final int nk = atomi.edges.size();
					if(nk>1) //spare the effort if it is a dead end
					{
						for(int k=0;k<nk;++k)
						{
							Bond edgek = atomi.edges.get(k);
						
							if( (edgek.left_ring == this) || (edgek.right_ring == this) )
							{	
								++shared_edges;
								/* //We actually dont relly on the correct order
								if( !has_vertex(atomi.neighbors.get(k)) )
								{
									System.out.println("ERROR while closing ring: " + rings.size() + " a bond is an edge but the Atoms are not both vertices");
								}
								*/ 
								if(!(has_vertex(edgek.a1) && has_vertex(edgek.a1)) )
								{
									System.out.println("Error while closing ring: " + rings.size() + " a bond is and edge but one of its atoms is no vertex");
								}
								if(!has_edge(edgek))
								{
									System.out.println("ERROR while closing ring: " + rings.size() + " a bond knowes the ring, but the ring does know the bond");
								}
							}
							else if(has_edge(edgek))
							{
								System.out.println("ERROR while closing ring: " + rings.size() + " the ring knows the bond, but the bond does know the ring");
							}
							
							if( ( !has_edge(edgek) ) && has_vertex(edgek.a1) && has_vertex(edgek.a2) )
							{
								++closing_edges;
							}
							
							if( (edgek.left_ring != null) && (edgek.right_ring != null) )
							{	++double_edges;}
							else if ( (edgek.left_ring != null)^(edgek.right_ring != null) )
							{   
								final int bridge_limit = Math.max(neighbors.size()-2,2);
								boolean bridge_head = (atomi.edges.size()<3);
								
								if(bridge_head)
								{
									int bridges = 0;
									int n = 0;
									while(n < nv)
									{	if(vertices.get(n++).edges.size()<3)
										{++bridges;}
									}
									bridge_head = (bridges < bridge_limit);
								}
								if(  bridge_head || 
								( (!neighbors.contains(edgek.left_ring) ) && (!neighbors.contains(edgek.right_ring)) ) )
								{	++single_edges;}
							}
							else
							{	++free_edges;}
						
						}
					}
					if( (shared_edges == 1) && ( (free_edges+single_edges > 0) || (closing_edges > 0) )  )
					{
						iend = i;
						break;	
					} 
				}
				if(iend == -1)
				{
					//supposed to happen at edges
					//System.out.println("could not find any candidate atom for growing ring " + rings.size() );
					return false;	
				}
				
				if(iend != 0 && iend != nv-1)
				{
					System.out.println("Warning, candidate atom for growing ring: " + rings.size() + " is not first or last vertex");
					debug_info();
				}
				 
				Atom atomj = vertices.get(iend);
				int kc = -1;
				double dsqr_min = Double.NaN;
				for(int k=0;k<atomj.edges.size();++k)
				{
					Atom atomk = atomj.neighbors.get(k);
					Bond bondk = atomj.edges.get(k);
					
					if(has_edge(bondk))
					{	continue;}
					
					boolean closing = ( has_vertex(bondk.a1) && has_vertex(bondk.a2) ); //(!has_edge(bondk)) &&
					
					if( (bondk.left_ring != null) && (bondk.right_ring != null) ) 
					{	
						if(debug_mode && closing)
						{
							System.out.println("Error cannot add closing edge to ring: " + rings.size() + 
							" because it is already winged by rings: " + bondk.left_ring.id + " & " + bondk.right_ring.id);
						}
						continue;
					} //ignore already known as well as fully winged edges
					
					
					Ring other = bondk.get_other_ring(null);
					
					if(other != null)
					{
						final double rrdistsqr = distsqr(pos,other.pos);
						final double rbdistsqr = distsqr(pos,bondk.pos); 
						if(rrdistsqr/rbdistsqr < 1.125 )
						{
							//System.out.println("rejecting bond behind a tentative neighbor, ring: " + rings.size() + " around: "  + (int)pos[0] + "," + (int)pos[1]);
							continue;
						}
					
					}
					 
					if( (!closing) && (other != null) && (atomj.neighbors.size()>2) && neighbors.contains(other))
					{	
						/*
						if(closing)
						{
							System.out.println("Error cannot add closing edge to ring: " + rings.size() + 
							" it would connect a second time to ring: " + other.id);
						}
						*/ 
						continue;
					} //dont walk along another ring if the vertex has three or more edges
					double dsqr = distsqr(pos,atomk.pos);
					if(closing)
					{
						dsqr = -1.0/(1.0+dsqr); //always prefer the shortest ring closing bond over all others 
					}
					// && (closing || ( atomk.edges.size() > atomk.rings.size() ) )
					if( ( !(dsqr >= dsqr_min) ) )
					{
						dsqr_min = dsqr;
						kc = k;
					}
				}
				if(kc != -1)
				{	
					if(!add_edge(atomj.edges.get(kc)))
					{
						return false;
					}
					
				}
				else
				{
					if(debug_mode)
					{
						System.out.println("could not pick any suitable bond of candidate atom for growing ring: " + rings.size() );
					}
					return false;
				}
				
			}
			/*
			if(!is_closed)
			{	
				System.out.println("Ups ring: " + rings.size() + " was supposed to be closed in line 1706");
				return false;
			}
			final int ns = neighbors.size();
			if(inverted && ns > 6)
			{
				System.out.println("Shit ring with " + ns + " neighbors reached line 1717");
			
				System.out.println("problematic mega ring id: " + rings.size() + " edges: " + edges.size() +
							 " vertices: " + vertices.size() + " neighbors: " + neighbors.size() +
							  " at x,y: " + (int)pos[0] + "," + (int)pos[1] + " radius: " + radius);
			
			}
			*/
			
			
			
			return true;
		
		}
		
		void make_springs()
		{
			final int ns = neighbors.size();
			for(int i = 0; i < ns; ++i)
			{
				final Ring ringi = neighbors.get(i);
				final Bond b = common_edge(ringi);
				if(b==null)
				{
					System.out.println("Error could not identify shared bond between rings: " + id + "&" + ringi.id );
					debug_info();
					ringi.debug_info();
					throw new RuntimeException("Bond common_edge(Ring other) failed");
				}
				else if( (b.spring == null) )
				{	new Spring(this,ringi);}
				
			}
			intensity = 0.0f;
			final int vs = vertices.size();
			for(int i = 0; i<vs; ++i)
			{	
				final Atom atomi = vertices.get(i);
				intensity -= atomi.intensity; //rings have opposite intensity		
			}
			intensity = intensity/vs*(0.5f*vs-1.0f); //average weight times full pies per polygon 
			
		}
		
		void sort_all()
		{
			//simple sort for theta of neighbors (rings)
			for(int i = 0; i<neighbors.size(); ++i)
			{	
				float[] veci = vecTo(pos,neighbors.get(i).pos);
				double thetai = Math.atan2(-veci[1],veci[0]);
				for(int j = i+1; j<neighbors.size(); ++j)
				{
					float[] vecj = vecTo(pos,neighbors.get(j).pos);
					final double thetaj = Math.atan2(-vecj[1],vecj[0]);
					if(thetaj < thetai)
					{
						Collections.swap(neighbors,i,j);
						thetai = thetaj;
					}
					
				}
			}
			//simple sort for theta of edges (bonds)
			for(int i = 0; i<edges.size(); ++i)
			{	
				float[] veci = vecTo(pos,edges.get(i).pos);
				double thetai = Math.atan2(-veci[1],veci[0]);
				for(int j = i+1; j<edges.size(); ++j)
				{
					float[] vecj = vecTo(pos,edges.get(j).pos);
					final double thetaj = Math.atan2(-vecj[1],vecj[0]);
					if(thetaj < thetai)
					{
						Collections.swap(edges,i,j);
						thetai = thetaj;
					}
					
				}
			}
			//simple sort for thetas of vertices (atoms)
			for(int i = 0; i<vertices.size(); ++i)
			{	
				float[] veci = vecTo(pos,vertices.get(i).pos);
				double thetai = Math.atan2(-veci[1],veci[0]);
				for(int j = i+1; j<vertices.size(); ++j)
				{
					float[] vecj = vecTo(pos,vertices.get(j).pos);
					final double thetaj = Math.atan2(-vecj[1],vecj[0]);
					if(thetaj < thetai)
					{
						Collections.swap(vertices,i,j);
						thetai = thetaj;
						
					}
				}
			}
		}
		
		private Ring(Bond b1, Bond b2)
		{
			if(b1 == null || b2 == null)
			{
				System.out.println("Cannot make a ring with null bonds");
				return; //dont register ourselves in rings
			}
			
			if(b1==b2)
			{
				System.out.println("caught double bond in new Ring");
				return; //dont register ourselves in rings
			}
			
			if(debug_mode)//this test does not scale to well with many rings
			{
				for(int i = 0; i < rings.size(); ++i)
				{
					Ring other = rings.get(i);
					if(other.has_edge(b1) && other.has_edge(b2))
					{
						//in principle that could become valid IF there are chains of atoms
						System.out.println("caught ambiguity (other ring has the same two bonds) in new Ring " + rings.size());
						System.out.println("offending ring");
						other.debug_info();
						System.out.println("this rings seeds: b1, b2");
						b1.debug_info();
						b2.debug_info();
						throw new RuntimeException("bad ring");
					}
				}
			} 
			if(b1.left_ring == null)
			{
				b1.left_ring = this;
			}
			else if(b1.right_ring == null)
			{
				b1.right_ring = this;
			}
			else
			{
				System.out.println("found already double used edge b1 in new Ring");
			}
			if( !edges.contains(b1))
			{	
				edges.add(b1);
				Ring other = b1.get_other_ring(this);
				
				if( (other != null) && (other != this) ) //second check is for vertices with only one bond 
				{
					if( !neighbors.contains(other) )
					{	
						neighbors.add(other);
						if(neighbors.size() > edges.size())
						{	System.out.println("too many neighbors (b1) in new Ring: " + rings.size());}
						
					}
					if( !other.neighbors.contains(this) )
					{	
						other.neighbors.add(this);
						if(other.neighbors.size() > other.edges.size())
						{	System.out.println("too many other.neighbors (b1) in new Ring: " + rings.size());}
					}
				}	
				
			}
			else
			{
				System.out.println("Error in new Ring() b1 was already known in edges");
			}
			
			
			
			if(b2.left_ring == null)
			{
				b2.left_ring = this;
			}
			else if(b2.right_ring == null)
			{
				b2.right_ring = this;
			}
			else
			{
				System.out.println("found already double used edge b2 in new Ring");
			}
			
			
			if( !edges.contains(b2))
			{	
				edges.add(b2);
				Ring other = b2.get_other_ring(this);
				if( (other != null) && (other != this) ) //this might happen if an vertex has only 2 edges
				{
					if( !neighbors.contains(other) )
					{	
						neighbors.add(other);	
						if(neighbors.size() > edges.size())
						{	System.out.println("too many neighbors (b2) in new Ring: " + rings.size());}	
						
					}
					if( !other.neighbors.contains(this) )
					{	
						other.neighbors.add(this);
						if(other.neighbors.size() > other.edges.size())
						{	System.out.println("too many other.neighbors (b2) in new Ring: " + rings.size());}
					}
				}
				
			}
			else
			{
				System.out.println("Error in new Ring() b2 was already known in edges");
			}
			
			
			
			if( !vertices.contains(b1.a1))
			{	
				vertices.add(b1.a1);
				b1.a1.rings.add(this);
				if( (b1.a1.edges.size() > 1) && (b1.a1.rings.size() > b1.a1.edges.size() ) )
				{
					System.out.println("Messed up topology of b1.a1 in new Ring: " + rings.size() );
				}	
			}
			if( !vertices.contains(b1.a2))
			{	
				vertices.add(b1.a2);
				b1.a2.rings.add(this);
				if((b1.a2.edges.size()>1) && (b1.a2.rings.size()>b1.a2.edges.size() ) )
				{
					System.out.println("Messed up topology of b1.a2 in new Ring: " + rings.size() );
				}		
			}
			if( !vertices.contains(b2.a1))
			{	
				vertices.add(b2.a1);
				b2.a1.rings.add(this);
				if((b2.a1.edges.size()>1) && (b2.a1.rings.size()>b2.a1.edges.size() ) )
				{
					System.out.println("Messed up topology of b2.a2 in new Ring: " + rings.size() );
				}	
			}
			if( !vertices.contains(b2.a2))
			{	
				vertices.add(b2.a2);
				b2.a2.rings.add(this);
				if( (b2.a2.edges.size()>1) && (b2.a2.rings.size()>b2.a2.edges.size() ) )
				{
					System.out.println("Messed up topology of b2.a2 in new Ring: " + rings.size() );
				}	
			}
			
			if(vertices.size() != 3)
			{
				System.out.println("found invalid number of inital vertices: " + vertices.size() + " in new Ring: " + rings.size() );
			}
			
			pos = Arrays.copyOf(vertices.get(0).pos,3);
			update_pos();
			/*
			//make_springs(); 
			//validate();
			if( (b1 == bonds.get(1586)) || (b2 == bonds.get(1586)) )
			{	System.out.println("new ring: " + rings.size() + " b1: " + b1.id + " b2: " + b2.id );}*/
			close();
			if( !is_closed )
			{
				
				if(debug_mode)
				{
					System.out.println("failed to close a new Ring: " + rings.size() + ", should only happen at edges");
					debug_info();
				}
				
				//make_springs(); 
				//validate();
				/*
				if(edges.contains(bonds.get(1586)))
				{debug_info(); }
				*/
				detach();
			}
			else
			{
				id = rings.size();
				rings.add(this);
				make_springs();
				validate();
				/*
				if(inverted && edges.size() > 4)
				{
					debug_info();
				}
				*/ 
				/*
				if(edges.contains(bonds.get(1586)))
				{debug_info(); }
				*/ 
				//System.out.println("established ring: " + id);
			}
		}
	
		void validate()
		{
			final int vs = vertices.size();
			for(int i = 0; i < vs; ++i)
			{
				boolean passed = true;
				Atom atomi = vertices.get(i);
				passed &= atomi.rings.contains(this);
				if(!passed) do_error("vertex: " + i );
			}
			
			
			final int es = edges.size();
			for(int i = 0; i < es; ++i)
			{
				boolean passed = true;
				Bond bondi = edges.get(i);
				passed &= ( (bondi.left_ring == this)^(bondi.right_ring == this) );
				
				
				
				passed &= ( (((bondi.left_ring == null)||(bondi.right_ring == null)) && (bondi.spring == null)) || 
							(((bondi.left_ring != null)&&(bondi.right_ring != null)) && (bondi.spring != null))		);
				//if( (bondi.spring != null) && (!springs.contains(bondi.spring)))
				if(!passed)
				{
					System.out.println("something is wrong with bond: " + bondi.id);
					bondi.debug_info();
					if(bondi.spring!=null)
					{	bondi.spring.debug_info();}
				}
				
				
				if(passed)
				{
					Ring other = bondi.get_other_ring(this);
					if(other != null)
					{
						passed &= neighbors.contains(other);
						if(bondi.spring != null)
						{
							passed &= springs.contains(bondi.spring);
							passed &= vertices.contains(bondi.spring.left_atom);
							passed &= vertices.contains(bondi.spring.right_atom);
							
							passed &= (bondi.spring.bond != null);
							passed &= edges.contains(bondi.spring.bond);
						
						}
						passed &= vertices.contains(bondi.a1);
						passed &= vertices.contains(bondi.a2);
					}
				}
				if(!passed) do_error("edge: " + i );
			}
			
			
			final int ns = neighbors.size();
			for(int i = 0; i < ns; ++i)
			{
				boolean passed = true;
				
				Bond bondi = common_edge(neighbors.get(i));
				passed &= (bondi != null);
				if(passed)
				{
					Spring springi = bondi.spring;
					passed &= springi != null;
					if(passed)
					{
						passed &= springs.contains(springi);
						passed &= springi.get_other_ring(this) == neighbors.get(i);
					}
				}
				if(!passed) do_error("neighbor: " + i );
			}
			
			
			final int ss = springs.size();
			for(int i = 0; i < ss; ++i)
			{
				boolean passed = true;
				Spring springi = springs.get(i);
				passed &= neighbors.contains(springi.get_other_ring(this));
				passed &= edges.contains(springi.bond);
				passed &= vertices.contains(springi.left_atom);
				passed &= vertices.contains(springi.right_atom);
				if(!passed) do_error("spring: " + i );
			}
		
			if(es > vs) do_error("edges: " + es + " > vertices:" + vs);
			if(ss != ns) do_error("springs: " + ss + " != neighbors:" + ns);
		
		}
	
		void do_error(String report)
		{
			System.out.println("Ring: " + (id==-1?rings.size():id) + " validation failed: " + report);
			debug_info();
			throw new RuntimeException(Macro.MACRO_CANCELED);
		}
	
	}
	
	public class Spring
	{
		public boolean detached = false;
		public float[] pos = null;
		public int id = -1;
		public Spring master = null;
		public int observers = 0;
		boolean is_interior = false;
		public Ring r1 = null;
		public Ring r2 = null;
		public Atom left_atom = null;
		public Atom right_atom = null;
		public Bond bond = null;
		private Line2D.Float[] segments = null;
		
		private Spring()
		{
			id = springs.size();
			springs.add(this);	
		}
		
		Atom get_other_atom(Atom a)
		{
			if(a == left_atom)
			{	return right_atom;}
			else if(a == right_atom)
			{	return left_atom;}
			else
			{	return null;}
			
		}
		
		
		public Ring get_other_ring(Ring r)
		{
			if(r == r1)
			{	return r2;}
			else if(r == r2)
			{	return r1;}
			else
			{	return null;}
			
		}
		
		float getLength()
		{
			return (float)Math.sqrt(distsqr(r1.pos,r2.pos));
		}
		
		boolean check_interior()
		{
			is_interior = (r1.is_interior) && (r2.is_interior);
			return is_interior;	
		}
		
		void addMark(Overlay overlay)
		{
			if(!detached)
			{
				Line2D.Float[] lines = get_lines(r1.pos,r2.pos,true);
				if(lines != null)
				{
					for(int i = 0; i < lines.length; ++i)
					{
						Line line = new Line ( lines[i].x1, lines[i].y1, lines[i].x2, lines[i].y2 );
						if(inverted)
						{
							line.setStrokeColor((is_interior||(!gray_edges))? Color.YELLOW.darker() : Color.GRAY.darker());
						}
						else
						{
							line.setStrokeColor((is_interior||(!gray_edges))? Color.GREEN.darker().darker() : Color.GRAY.darker());
						}
						if(impSlice != 0)
						{	line.setPosition(impSlice);}
						overlay.add( line,"Spring:" + id);
					}
				}							
			}
		}
		
		public boolean intersecting() 	
		{	
			final int ss = springs.size();
			final Line2D.Float[] msegs = segments;
			for(int i = 0; (i < ss); ++i)
			{	
				Spring spi = springs.get(i);
				if(spi.r1 == r1 || spi.r2 == r2 || spi.r1 == r2 || spi.r2 == r1 )
				{	continue;}
				final Line2D.Float[] osegs = spi.segments;	
				for(int m = 0; m < msegs.length; ++m)
				{
					final Line2D.Float ms = msegs[m];
					for(int o = 0; o < osegs.length; ++o)
					{
						if( ms.intersectsLine(osegs[o])) return true;
					}
				}
			}	
			return false;
		}
		
	
		private Spring(Ring rr1, Ring rr2)
		{
			r1 = rr1;
			r2 = rr2;
			segments = get_lines(r1.pos,r2.pos,false);
			//if(!intersecting()) //this always works out if the rings are meaningful
			boolean first = true;
			final int es = r1.edges.size();
			for(int i=0; i < es; ++i)	
			{	
				final Bond bondi = r1.edges.get(i);
				
				if(r2.has_edge(bondi))
				{
					bondi.spring = this;
					if(first)
					{
						bond = bondi;
						left_atom = bondi.a1;
						right_atom = bondi.a2;
					}
				}
				
			}	
			id = springs.size();
			springs.add(this);
			
			r1.springs.add(this);
			r2.springs.add(this);
				
		}
	
		void debug_info()
		{
			System.out.println("spring: " + id + " in mesh: " + impSlice  +" detached: " + detached);
			System.out.println("rings: " + r1.id + "&" + r2.id );
			System.out.println("bond: " + bond.id + " left&right atoms: " + left_atom.id + "&" + right_atom.id);	
		}
	
	
	
	}
	
	
	public class Signature
	{
		public Atom atom = null;
		public String sig_string = null;
		
		private Signature(String master_sig_string)
		{
			if(hex_pixels)
			{	throw new RuntimeException("You cannot use signatures with hexagonal pixels");}
			sig_string = master_sig_string;
		}
		
		
		private Signature(Atom a0)
		{
			if(hex_pixels)
			{	throw new RuntimeException("You cannot use signatures with hexagonal pixels");}
			this.atom = a0;
			
			
			
			
			int p5 = -1; //position of 5-ring
			int p6 = -1;
			int p7 = -1; //position of 7-ring
			final int ars = atom.rings.size();
			for(int i = 0; i < ars; ++i)
			{
				 final int ns = atom.rings.get(i).neighbors.size();
				 switch(ns)
				 {
					 case 5: p5=i; break;
					 case 6: p6=i; break;
					 case 7: p7=i; break;
					 default:
				 }
			}
			if((ars!=3) || (p5==-1) || (p6==-1) || (p7==-1) )
			{
				sig_string = "unspecific";
				return;
			}
			
			ArrayList<Ring> known_rings = new ArrayList<Ring>(rings.size()/4);
			ArrayList<Ring> fringe = new ArrayList<Ring>(rings.size()/4);	
			StringBuilder sig_builder = new StringBuilder(256);
			boolean fringe_is_interior = false;
			
			
			fringe.add(atom.rings.get(p5));
			fringe.add(atom.rings.get(p6));
			fringe.add(atom.rings.get(p7));
			
			
			//HAS DUPLICATE a bottom of loop
			int fs = fringe.size();
			fringe_is_interior = ( (fs==0) ? false : fringe.get(0).is_interior);
			for(int i = 1; (i < fs && fringe_is_interior); ++i)
			{
				fringe_is_interior &= ( fringe.get(i).is_interior );
			}
			int lvl = 1;
			while(fringe_is_interior)
			{
				//add a line to string representation
				for(int i = 0; (i < fs); ++i)
				{
					final Ring fringi = fringe.get(i);
					int n = fringe.get(i).edges.size();
					sig_builder.append(n);
					if( (i+1) < fs )
					{	sig_builder.append('-');}
					//System.out.print("add " + n + "-ring at: " + (int)fringi.pos[0] + "," + (int)fringi.pos[1] + "\t");
				}
				sig_builder.append('\n');
				sig_string = sig_builder.toString(); //store and save what we have safely identified so far
				//System.out.println("\ncompleted fringe: " + lvl + " rings: " + fs);
				ArrayList<Ring> nfringe = new ArrayList<Ring>(rings.size()/4);
				
				//process the recent fringe
				Ring fringin = fringe.get(fs-1);
				for(int i = 0; (i < fs); ++i)
				{
					final Ring fringi = fringe.get(i);
					//walk around ringi starting from either initial fringin or the last fringin
					do
					{
						final Bond brad = get_Bond(fringi,fringin); //a radial bond reaching the outer edge of the fringe 
						Atom a1 = null;
						Atom a2 = null;
						
						if(brad == null) //that will be necessary if there are 4-fold coordinated atoms
						{
							for(int m = 0; m < fringi.vertices.size(); ++m)
							{
								final Atom am = fringi.vertices.get(m);
								for(int n = 0; n < fringin.vertices.size(); ++n)
								{
									final Atom an = fringin.vertices.get(n);
									if( (am == an) && (am.rings.size()>2) )
									{
										if(a1 == null)
										{	a1 = am;}
										else if(a2 == null)
										{	a2 = am;}
										else
										{	
											System.out.println("serious problem two rings share more than two atoms with at least 3 neighbors");
											fringi.debug_info();
											fringin.debug_info();
											throw new RuntimeException(Macro.MACRO_CANCELED);	
										}
									}
								}
							}	
						}
						else
						{
							a1 = brad.a1;
							a2 = brad.a2;
						}
						if(a1 == null) //that should only happen at the edges
						{
							/*
							System.out.println("serious problem could not find a a shared atom between");
							fringi.debug_info();
							fringin.debug_info();
							throw new RuntimeException(Macro.MACRO_CANCELED);
							*/
							
							return;
							
							 
						}
						
						Ring rm = null;
						//look for a yet unknown ring that shares an edge with fringin and an atom with fringi 
						int maxk = a1.rings.size();
						if(a2!=null && a2.rings.size() > maxk)
						{	maxk = a2.rings.size();}
						for(int k = 0; k < maxk;++k)
						{
							Ring rk = null;
							if( (a1!=null) && (k < a1.rings.size()) )
							{
								rk = a1.rings.get(k);
								if( (rk != fringi) && (rk != fringin) && fringin.neighbors.contains(rk) && 
									(!fringe.contains(rk)) && (!nfringe.contains(rk)) && (!known_rings.contains(rk)))
								{
									rm = rk; 
									break;	
								}
							}
							if( (a2!=null) && (k < a2.rings.size()) )
							{
								rk = a2.rings.get(k);
								if( (rk != fringi) && (rk != fringin) && fringin.neighbors.contains(rk) && 
									(!fringe.contains(rk)) && (!nfringe.contains(rk)) && (!known_rings.contains(rk)))
								{
									rm = rk; 
									break;	
								}
							}		
						}
						if(rm == null) {break;}// we walked the open perimeter
						nfringe.add(rm);
						fringin = rm;
					}while(true);
					/*
					final int ns = fringi.neighbors.size();
					
					//determine the index of the first already known neighboring ring 
					int j0 = 0;
					for(int j = 0; j < ns; ++j)
					{
						Ring rj = fringi.neighbors.get(j);
						if( fringe.contains(rj) || nfringe.contains(rj) || known_rings.contains(rj) )
						{	j0 = j;	break;}
					}
					
					
					
					
					
					//put neigboring rings of frings to nfring in counterclockwise order
					//as long as their clockwise predecesor is already known
					for(int j = 0; j < ns; ++j)
					{
						final int n = (j+j0)%ns;
						final int m = (n+1)%ns;
						final Ring rn = fringi.neighbors.get(n); 
						final Ring rm = fringi.neighbors.get(m);	
						if( (fringe.contains(rn) || nfringe.contains(rn) || known_rings.contains(rn)) && 
							(!fringe.contains(rm)) && (!nfringe.contains(rm)) && (!known_rings.contains(rm) ) )
						{	nfringe.add(rm);}
					}
					*/	
				}
				
				//version for softies, but actually needed with periodic meshes
				//put the processed fringe to the core of known_rings
				/*
				for(int i = 1; (i < fs); ++i)
				{
					Ring ringi = fringe.get(i);
					known_rings.add(ringi);	
				}
				*/
				known_rings = fringe;
				
				fringe = nfringe;
				//DUPLICATE with loop initialization
				fs = fringe.size();
				fringe_is_interior = ( (fs==0) ? false : fringe.get(0).is_interior);
				for(int i = 1; (i < fs && fringe_is_interior); ++i)
				{
					fringe_is_interior &= ( fringe.get(i).is_interior );
				}
				++lvl;
				
			}
		}
		
		//give positive score for bigger and bigger matching neighborhood
		//-1 for a confirmed deviation on Rings that are interior 
		public int matchTo(Atom oa)
		{
			if(rate_signature(atom) != 2)
			{	return -1;}
			Signature other = new Signature(oa);
			return eval_sig_strings(sig_string,other.sig_string);
		}
	}
	
	Bond get_Bond(Atom a1, Atom a2) //local search
	{
		for(int i=0; i < a1.neighbors.size(); ++i)
		{
			if(a1.neighbors.get(i)==a2)
			{	return a1.edges.get(i); }		   	    
		}
		return null;
	}
	
	Bond get_Bond(Ring r1, Ring r2) //search for the one shared bond
	{
		Bond bu = null;
		for(int i=0; i < r1.edges.size(); ++i)
		{
			Bond b = r1.edges.get(i);
			if(r2.has_edge(b) )
			{	
				if(bu == null)
				{
					bu = b;
				}
				else
				{
					bu = null;
					break;
				}
			}    
		}
		return bu;
	}
	
	
	
	Atom get_Atom(Atom a1, Atom a2) //local search 
	{
		for(int i = 0; i < a1.neighbors.size(); ++i)
		{
			if(a2.neighbors.contains(a1.neighbors.get(i)))
			{
				return a1.neighbors.get(i);
			}
		}
		return null;
	}
	
	
	
	Ring get_Ring(Bond b1, Bond b2) //local search
	{
		if( (b1.left_ring!=null) &&
			( (b1.left_ring == b2.left_ring) || 
			  (b1.left_ring == b2.right_ring) ) 
		  )
		{	return b1.left_ring;}
		else if( (b1.right_ring!=null) &&
			( (b1.right_ring == b2.left_ring) || 
			  (b1.right_ring == b2.right_ring) ) 
		  )
		{	return b1.right_ring;}
		
		return null;
	}
	
	public double distsqr(float[] p1, float[] p2)
	{
		float dx = p2[0]-p1[0];
		float dy = p2[1]-p1[1];
		float dz = p2[2]-p1[2];
		if(!hex_pixels)
		{	
			if( periodic)
			{
				dx = (dx+3*impWidth/2)%impWidth-impWidth/2;
				dy = (dy+3*impHeight/2)%impHeight-impHeight/2;
			}
			return( dx*dx + dy*dy + dz*dz );
		}
		else
		{
			int[] xyz = new int[]{Math.round(dx),0,Math.round(dz)};
			xyz[1] = -xyz[0] -xyz[2];
			imgHP.periodic_xyz(xyz);
			return (double)(xyz[0]*xyz[0] + xyz[2]*xyz[2] + xyz[0]*xyz[2]); 
		}
	}
	
	public float[] vecTo(float[] p1, float[] p2)
	{
		float dx = p2[0]-p1[0];
		float dy = p2[1]-p1[1];
		float dz = p2[2]-p1[2];
		if(!hex_pixels)
		{
			if(periodic)
			{
				dx = (dx+3*impWidth/2)%impWidth-impWidth/2;
				dy = (dy+3*impHeight/2)%impHeight-impHeight/2;
			}
			
			float[] v12 = {dx,dy,dz};
			return v12;
		}
		else
		{
			int[] xyz = new int[]{Math.round(dx),0,Math.round(dz)};
			xyz[1] = -xyz[0] -xyz[2];
			imgHP.periodic_xyz(xyz);
			return new float[]{ (float)xyz[0], (float)xyz[1], (float)xyz[2] };
		}
	}
	
	public double dist2Dsqr(float[] p1, float[] p2)
	{
		if(!hex_pixels)
		{
			float dx = p2[0]-p1[0];
			float dy = p2[1]-p1[1];
			if( periodic)
			{
				dx = (dx+3*impWidth/2)%impWidth-impWidth/2;
				dy = (dy+3*impHeight/2)%impHeight-impHeight/2;
			}
			return (double)( dx*dx + dy*dy );
		}
		else //hex_pixels are always in cube coordinates anyways
		{	return distsqr(p1,p2);} 
	}
	
	public Atom put_atom(float[] pos)
	{
		if(hex_pixels)
		{
			int[] cube = ixyz(pos);
			imgHP.periodic_xyz(cube);
			pos[0] = cube[0];
			pos[1] = cube[1];
			pos[2] = cube[2];
		}
		else if(periodic)
		{
			pos[0] = (pos[0]+impWidth)%impWidth;
			pos[1] = (pos[1]+impHeight)%impHeight;
		}
		
		return new Atom(pos);
	}
	
	private Atom put_atom()
	{
		return new Atom();
	}
	
	public void remove_atom(Atom a)
	{	
		a.detach();
	}
	
	private Bond put_bond()
	{
		return new Bond();
	}
	
	private Ring put_ring()
	{
		return new Ring();
	}
	
	public void refresh()
	{
		
		if(!(needs_update))
		{	
			return;
		}
		removefrom_master(); // abandon all master relations, the new topology might differ
		
		ArrayList<Atom> atoms2 = new ArrayList<Atom>(atoms.size());
		for(int i = 0; i < atoms.size(); ++i)
		{
			final Atom atomi = atoms.get(i);
			if(!atomi.detached)
			{	
				atomi.id = atoms2.size();
				atoms2.add(atomi);	
			}
		}
		atoms = atoms2; //purge all detached atoms
		fatoms = atoms.toArray( new Atom[0] );
		bonds.clear();
		fbonds = null;
		rings.clear();
		frings = null;
		springs.clear();
		fsprings = null;
		
		
		int candidates = 0;
		for(int i = 0; i < fatoms.length; ++i )
		{
			Atom atomi = fatoms[i];
			atomi.rings.clear();
			atomi.edges.clear();
			atomi.neighbors.clear();
			atomi.candidates.clear();
			atomi.shortest_distsqr = Float.NaN;
			atomi.short_distsqr = Float.NaN;
			atomi.long_distsqr = Float.NaN;
			for(int j = 0; j < i; ++j )
			{
				Atom atomj = fatoms[j];
				atomi.addCandidate(atomj);
			}
		}
		int num_bonds = bonds.size();
		do
		{	
			num_bonds = bonds.size();
			for(int i = 0; i < fatoms.length; ++i )
			{	
				Atom atomi = fatoms[i];
				final int atc = atomi.top_candidate;
				Atom atomj = (atc!=-1)?atomi.candidates.get(atc):null;
				final int btc = (atc == -1) ? -1 : atomj.top_candidate;
				if( (atc != -1) && (btc != -1) &&
					 atomi.short_distsqr <= atomi.distsqrTo(atomj)) //*(inverted?1.1:1.0)
				{
					atomi.addNeighbor(atomj); 
				}
			}
		
		}while(bonds.size() > num_bonds);
		
		ArrayList<Bond> bonds2 = new ArrayList<Bond>(bonds.size());
		for(int i = 0; i < bonds.size(); ++i)
		{
			final Bond bondi = bonds.get(i);
			if(!bondi.detached)
			{	
				bondi.id = bonds2.size();
				bonds2.add(bondi);	
			}
		}
		bonds = bonds2; //purge all detached bonds
		
		
		
		fbonds = bonds.toArray( new Bond[0] );
		
		
		for(int i=0; i<fatoms.length; ++i)
		{	
			fatoms[i].sort_neighbors();
			fatoms[i].check_interior();	
		}
		interior_bonds = 0;
		for(int i = 0; i < fbonds.length; ++i)
		{	
			interior_bonds += (fbonds[i].check_interior()?1:0);
		}
		
		
		num_unique_pairs = 0;
		for(int i = 0; i < fatoms.length; ++i)
		{	num_unique_pairs += fatoms[i].init_unique_pairs();}
		
		boolean[] trial_mask = new boolean[fatoms.length]; // all false
		int ti = -1;
		int max_pick_errors = 5;
		int num_rings = 0;
		
		do
		{
			num_rings = rings.size();
			ti = -1;
			int tightest_spot = 0;
			int highest_coordination = 0;
			final float[] cpos = hex_pixels ? new float[]{0.0f,0.0f,0.0f} : new float[]{0.5f*impWidth,0.5f*impHeight, 0.0f};
			double min_distsqr = Double.NaN;
			for(int i = 0; i < fatoms.length; ++i)
			{
				if(trial_mask[i] )
				{	continue;}
				final int es = fatoms[i].edges.size();
				if(es < 3)
				{	continue;}
				final int rs = fatoms[i].rings.size();
				if( (es>1) && (rs > es) )
				{
					System.out.println("bad atom near x,y: " + (int)fatoms[i].pos[0] + "," + (int)fatoms[i].pos[1] +
					" neighbors: " + fatoms[i].neighbors.size() + " edges: " + fatoms[i].edges.size() + 
					" rings: " + fatoms[i].rings.size());
					throw new RuntimeException("encountered more rings than edges at a vertex");
				}
				final double dsqr = dist2Dsqr(cpos,fatoms[i].pos);
				//select by a highest coordination
				//			b least missing rings
				//			c closest to center
				//          The priority matters as Ring close cannot properly handle some unlucky seeds at the fringe
				
				if( (es == rs) || (es<2) )
				{	continue;}
				boolean take = false;
				boolean on_tie = false;
				take = ( (es>2) && (es>highest_coordination) );
				on_tie = (es == highest_coordination);
				if(!(on_tie || take)) continue;
				
				if(!take)
				{
					take = ( (es-rs > tightest_spot) );
					on_tie = ( (es-rs) == tightest_spot);
					if(on_tie)
					{
						take = (!(dsqr > min_distsqr) );
					}
				}
				
				if( take )
				{
					ti=i;
					tightest_spot = es-rs;
					min_distsqr = dsqr;
					highest_coordination = es;
				}
			}
			
			if(ti != -1)
			{		
				Atom atomi = fatoms[ti];
				trial_mask[ti] = true;
				int i1 = -1;
				int i2 = -1;
				final int es = atomi.edges.size();
				int trials_per_atom = ((es>2)?es:0);
				do
				{		
					if(--trials_per_atom < 0)
					{	break;}
					i1 = -1;
					i2 = -1; 
					//look for two consequtive not fully winged bonds
					for(int i0 = 0; i0 < es; ++i0)
					{
						Bond b1 = atomi.edges.get(i0);
						Bond b2 = atomi.edges.get((i0+1)%es);
						if( ( (b1.left_ring==null) || ( b1.right_ring == null) ) &&
							( (b2.left_ring==null) || ( b2.right_ring == null) ) )		
						{
							Ring b1r = ( (b1.left_ring==null) ? b1.right_ring : b1.left_ring );
							Ring b2r = ( (b2.left_ring==null) ? b2.right_ring : b2.left_ring );
							if( ( (b1r == null) || (b2r == null) ) || (b1r != b2r) )
							{
								i1 = i0;
								i2 = (i0+1)%es;
								//System.out.println("picked atom: " + atomi.id + " ,11: " + i1 + " ,12: " + i2 + " for next ring: " + rings.size());
								break;
							}
						}
					}
					/*
					if( !( (i1 != i2) && (i1!=-1) && (i2!=-1) ) ) //hmm just fall back to possibly unordered bonds/neighbors
					{
						System.out.println("prescan found no more ordered starting bonds for ring: " + rings.size() + " near atom" );
						atomi.debug_info();
						//pick fully unwinged edges
						for(int i0 = 0; i0 < es; ++i0)
						{
							Bond bondi0 = atomi.edges.get(i0);
							if( (bondi0.left_ring==null) && ( bondi0.right_ring == null) )
							{
								if(i1 == -1)
								{	i1 = i0;}
								else if(i2 == -1)
								{ i2 = i0;}
							}
						}
						//finally fill up with single winged edges
						for(int i0 = 0; i0 < es; ++i0)
						{
							Bond bondi0 = atomi.edges.get(i0);
							if( (bondi0.left_ring==null)^( bondi0.right_ring == null) )
							{
								if(i1 == -1)
								{	i1 = i0;}
								else if(i2 == -1)
								{ i2 = i0;}
							}
						}
						if( (i1 != -1) && (i2 != -1) )
						{
							System.out.println("picking unordered bonds:");
							atomi.edges.get(i1).debug_info();
							atomi.edges.get(i1).debug_info();
						}	
					}
					*/
					if( (i1 != i2) && (i1!=-1) && (i2!=-1) )
					{
						
						if( ((i1+1)%atomi.edges.size() != i2) && (atomi.edges.size()>3))
						{
							//that should be taken care of by Ring.cut_tail()
							System.out.println("Attention: ring: "+ rings.size()  +" will be initialized with non consequtive edges of atom: " + atomi.id);
							atomi.debug_info();	
						}
						
						Bond b1 = atomi.edges.get(i1);
						Bond b2 = atomi.edges.get(i2);
						Ring b1r = ( (b1.left_ring==null) ? b1.right_ring : b1.left_ring );
						Ring b2r = ( (b2.left_ring==null) ? b2.right_ring : b2.left_ring );
						if( ( (b1r == null) || (b2r == null) ) || (b1r != b2r) )
						{
							//System.out.println("seeding ring: " +rings.size());
							new Ring(b1,b2);
						}
						else
						{
							System.out.println("failed to predict any suitable bonds for new ring: " + rings.size());
							ti = -1;
						}	
					}
				}
				while( (i1!=-1) && (i2!=-1) && (i1!=i2) );		
			}	
		}
		while( (ti != -1) ); 
		
		
		ArrayList<Ring> rings2 = new ArrayList<Ring>(rings.size());
		for(int i = 0; i < rings.size(); ++i)
		{
			final Ring ringi = rings.get(i);
			if(!ringi.detached)
			{	
				ringi.id = rings2.size();
				rings2.add(ringi);
				if(inverted && ringi.springs.size() > 10)
				{
					System.out.println("Error: found surviving monster ring:");
					ringi.debug_info();
				}
					
			}
		}
		rings = rings2; //purge all detached rings
		
		
		
		frings = rings.toArray( new Ring[0] );
		for(int i = 0; i < frings.length; ++i)
		{
			frings[i].check_interior();
			frings[i].sort_all();	
		}
		fsprings=springs.toArray( new Spring[0]);
		for(int i = 0; i < fsprings.length; ++i)
		{
			fsprings[i].check_interior();
		}
		needs_update = false;
	}

	//simply accept any external topological manipulations
	public void fake_refresh() 
	{
		
		
		ArrayList<Atom> atoms2 = new ArrayList<Atom>(atoms.size());
		for(int i = 0; i < atoms.size(); ++i)
		{
			final Atom atomi = atoms.get(i);
			if(!atomi.detached)
			{	
				atomi.id = atoms2.size();
				atoms2.add(atomi);	
			}
		}
		atoms = atoms2; //purge all detached atoms
		fatoms = atoms.toArray( new Atom[0] );
		
		ArrayList<Bond> bonds2 = new ArrayList<Bond>(bonds.size());
		for(int i = 0; i < bonds.size(); ++i)
		{
			final Bond bondi = bonds.get(i);
			if(!bondi.detached)
			{	
				bondi.id = bonds2.size();
				bonds2.add(bondi);	
			}
		}
		bonds = bonds2; //purge all detached bonds
		
		
		
		fbonds = bonds.toArray( new Bond[0] );
		
		ArrayList<Ring> rings2 = new ArrayList<Ring>(rings.size());
		for(int i = 0; i < rings.size(); ++i)
		{
			final Ring ringi = rings.get(i);
			if(!ringi.detached)
			{	
				ringi.id = rings2.size();
				rings2.add(ringi);		
			}
		}
		rings = rings2; //purge all detached rings
		frings = rings.toArray( new Ring[0] );
		
		
		ArrayList<Spring> springs2 = new ArrayList<Spring>(springs.size());
		for(int i = 0; i < springs.size(); ++i)
		{
			final Spring springi = springs.get(i);
			if(!springi.detached)
			{	
				springi.id = springs2.size();
				springs2.add(springi);		
			}
		}
		springs = springs2; //purge all detached springs
		fsprings = springs.toArray( new Spring[0] );
		
		needs_update = false;
	}


	public void drawMarks(ImagePlus bild, boolean draw_atoms, boolean draw_bonds, boolean draw_rings, boolean draw_springs)
	{
		Overlay overlay = new Overlay();
		if(!needs_update || (this == master_mesh) )
		{
			if(draw_atoms)
			{
				for(int i=0; i < fatoms.length; ++i)
				{	fatoms[i].addMark( overlay );}
			}
			
			if(draw_bonds)
			{
				for(int i=0; i < fbonds.length; ++i)
				{	fbonds[i].addMark( overlay );}
			}
			
			if(draw_rings)
			{
				for(int i=0; i < frings.length; ++i)
				{	frings[i].addMark( overlay );}
			}
			
			if(draw_springs)
			{
				for(int i=0; i < fsprings.length; ++i)
				{	fsprings[i].addMark( overlay );}
			}
			 
		}
		overlay.translate(0.5,0.5);
		if(this == master_mesh)
		{
			bild.setOverlay(overlay);
			return;
		}
		
		
		Overlay previous = bild.getOverlay();
		if(previous != null)
		{
			
			Roi[] prois = previous.toArray();
			Roi[] mrois = overlay.toArray();
			Roi[] jrois = new Roi[prois.length + mrois.length];
			
			int s = 0;
			for(int i = 0; i < prois.length; ++i)
			{
				if(prois[i]!=null && prois[i].getPosition()!=impSlice)
				{	jrois[s++] = prois[i];}
			}
			for(int i = 0; i < mrois.length; ++i)
			{	
				if(mrois[i]!=null)
				{	jrois[s++] = mrois[i];}
			}
			Overlay jol = new Overlay();
			for (int k = 0; k < s; ++k)
			{	jol.add(jrois[k]);}
			
			bild.setOverlay(jol);
		}
		else
		{	bild.setOverlay(overlay);}
	}


	public Signature create_signature()
	{
		System.out.println("creating signature " + impSlice + "...");
		final float[] cpos = hex_pixels ? new float[]{0.0f,0.0f,0.0f} : new float[]{0.5f*impWidth,0.5f*impHeight, 0.0f};
		
		Atom best_atom = null;
		double dsqr_min = Double.NaN;
		for(int i = 0; i < fatoms.length; ++i)
		{
			Atom atomi = fatoms[i];
			final double dsqr = distsqr(cpos,atomi.pos);
			if( (!(dsqr_min < dsqr)) && (rate_signature(atomi)>1) )
			{
				dsqr_min = dsqr;
				best_atom = atomi;
			}
		}
		if(best_atom == null)
		{
			System.out.println("failed to pick any suitable 5-6-7 atom in mesh: " + impSlice);
                        return null;
		}
		System.out.println("picked atom at x,y: " + (int)best_atom.pos[0] + "," + (int)best_atom.pos[1]);
		signature = new Signature(best_atom);
		//System.out.print(signature.sig_string);
		return signature;
	}
	
	public Signature link_signature(String sig_str1)
	{
		if(this == master_mesh && atoms.isEmpty())
		{
			signature = new Signature(sig_str1);
			return signature;
		}
		
		//System.out.println("linking atoms of slice " + impSlice  + " to signature ...");
		Signature best_sig = null;
		int max_score = 0;
		for(int i = 0; i < fatoms.length; ++i)
		{
			Atom atomi = fatoms[i];
			if((rate_signature(atomi)>1) )
			{
				Signature tmp_sig = new Signature(atomi);
				final int score = eval_sig_strings(sig_str1, tmp_sig.sig_string);
				System.out.println("checked 5-6-7 atom at: " +(int)atomi.pos[0] + "," + (int)atomi.pos[1] + " score: " + score);
				if(score > max_score)
				{
					max_score = score;
					best_sig = tmp_sig;
				}
				else if(score == max_score)
				{
					if(IJ.showMessageWithCancel("Ambiguity", "Is the latest point the best anchor?"))
					{
						best_sig = tmp_sig;
						++max_score;
					}
				}
			}
		}
		if(best_sig == null)
		{
			System.out.println("failed to pick matching 5-6-7 atom in mesh: " + impSlice);
			return null;
		}
		signature = best_sig;
		//System.out.println("picked atom at x,y: " + (int)signature.atom.pos[0] + "," + (int)signature.atom.pos[1]);
		//System.out.print(signature.sig_string);
		return signature;	
	}
	
	public void incorporate(boolean add) //false means only removing
	{
		if(atoms.isEmpty())
		{	return;}
		
		if(add)
		{
			if(master_mesh.signature != null)
			{	signature = link_signature(master_mesh.signature.sig_string);}
			else
			{
				if(create_signature()==null){ return;}
				master_mesh.link_signature(signature.sig_string);	
			}
			if(signature == null)
			{	return;}
			addto_master();
			clean_up_master_mesh();
		}
		else
		{
			removefrom_master();
		}
		
		if(master_mesh.atoms.size()>12) //actually many more than 1 are required
		{	master_mesh.create_signature();}
		
	}
	
	private void addto_master()
	{
		
		Atom sa = signature.atom;
		Atom sam = sa.master;
		Atom msa = master_mesh.signature.atom;
		if( (sam != null) && (msa != null) && (sam != msa) )
		{
			System.out.println("Warning surfacemesh " + impSlice + " has changed primary matching atom to master_mesh");
			removefrom_master();
			//clean_up_master_mesh();
			sam = null;
		}
		if( (sam == null) && (msa != null) )
		{
			sam = msa; //the actual sa.master is NOT affected by that
			
			int mi5 = -1;
			int mi6 = -1;
			int mi7 = -1;
			int i5 = -1;
			int i6 = -1;
			int i7 = -1;
			
			for(int i = 0; i < 3; ++i) 
			{
				switch( msa.rings.get(i).vertices.size() ) 
				{
					case 5: mi5 = i; break;
					case 6: mi6 = i; break;
					case 7: mi7 = i; break;
					default:
				}
				switch( sa.rings.get(i).vertices.size() ) 
				{
					case 5: i5 = i; break;
					case 6: i6 = i; break;
					case 7: i7 = i; break;
					default:
				}
			}
			Ring r5 = sa.rings.get(i5);
			Ring mr5 = msa.rings.get(mi5);
			Ring r7 = sa.rings.get(i7);
			Ring mr7 = msa.rings.get(mi7);
			Bond b57 = get_Bond(r5,r7);
			Bond mb57 = master_mesh.get_Bond(mr5,mr7);
			
			r5.master = mr5;
			r5.master.observers++;
			b57.master = mb57;
			b57.master.observers++;
			sa.master = msa;
			sa.master.observers++;
			
			Bond nb = sa.get_other_bond(r5,b57);
			Bond mnb = msa.get_other_bond(mr5,mb57);
			Atom na = nb.get_other_atom(sa);
			Atom mna = mnb.get_other_atom(msa);
			//TODO check that none of them is null
			int num_bonds = 1;
			int num_atoms = 1;
			while( (nb.master == null) || (na.master == null) )
			{
				if(nb.master == null)
				{
					nb.master = mnb;
					nb.master.observers++;
					++num_bonds;
				}
				else if(nb.master != mnb)
				{
					System.out.println("Shit, nb.master != mnb while placing central 5-ring");
				}
				if(na.master == null)
				{
					na.master = mna;
					na.master.observers++;
					++num_atoms;
				}
				else if(na.master != mna)
				{
					System.out.println("Shit, na.master != mna while placing central 5-ring");
				}
				
				nb = na.get_other_bond(r5,nb);
				mnb = mna.get_other_bond(mr5,mnb);
				na = nb.get_other_atom(na);
				mna = mnb.get_other_atom(mna);	
			} 
			if( (num_bonds != 5) || (num_atoms != 5) )
			{
				System.out.println("Ups, just put a 5-ring and got " + num_atoms + " vertices and " + num_bonds + " edges" );
			}
			//push_atom_relations();
			//r5.master.sort_all();
			
			if( !validate_master_ring(r5) )
			{
				System.out.println("Ups something went wrong identifieying the first ring on an already populated joined mesh");
				r5.debug_info();
				r5.master.debug_info();
				removefrom_master(); 
				//clean_up_master_mesh();
				return;
			}
	
		}
		
		
		if(sam == null) //provide a seeding ring that will contain the signature.atom anyways
		{ //no checks here as it ought to be part of a valid signature anyways
			if(master_mesh.atoms.size() > 0)
			{
				System.out.println("Error mastermesh containes atoms, but signature atoms did not lign up");
			}
			
			Ring r0 = sa.rings.get(0); 
			r0.master = master_mesh.put_ring();
			r0.master.observers = 1;
			r0.master.pos = Arrays.copyOf(r0.pos,3);
			r0.master.is_interior = true;
			r0.master.is_closed = true;
			r0.master.radius = r0.radius;
			for(int i = 0; i < r0.vertices.size(); ++i)
			{
				Atom ai = r0.vertices.get(i);
				ai.master = master_mesh.put_atom();
				ai.master.observers = 1;
				ai.master.pos = Arrays.copyOf(ai.pos,3);
				ai.master.is_interior = ai.is_interior;
				ai.master.shortest_distsqr = ai.shortest_distsqr;
				r0.master.vertices.add(ai.master);
				ai.master.rings.add(r0.master);
			}
			for(int i = 0; i < r0.edges.size(); ++i)
			{
				Bond bi = r0.edges.get(i);
				bi.master = master_mesh.put_bond();
				bi.master.observers = 1;
				bi.master.pos = Arrays.copyOf(bi.pos,3);
				bi.master.is_interior = bi.is_interior;
				r0.master.edges.add(bi.master);
				bi.master.left_ring = r0.master;
				bi.master.a1 = bi.a1.master;
				bi.master.a2 = bi.a2.master;
				
				bi.master.a1.neighbors.add(bi.master.a2);
				bi.master.a2.neighbors.add(bi.master.a1);
				bi.master.a1.edges.add(bi.master);
				bi.master.a2.edges.add(bi.master);	
			}
			//push_atom_relations();
			r0.master.sort_all();
			sam = sa.master;
			//r0.master still has no neighbors	
			
			if( !validate_master_ring(r0) )
			{
				System.out.println("Ups created a faulty seed ring - master pair");
				r0.debug_info();
				r0.master.debug_info();
				removefrom_master(); 
				//clean_up_master_mesh();
				return;
			}
		
		
		}
		
		
		float[] shift = new float[3];
		shift[0] = sam.pos[0]-sa.pos[0];
		shift[1] = sam.pos[1]-sa.pos[1];
		shift[2] = sam.pos[2]-sa.pos[2];
                
                manual_init();
		/*
		System.out.println("shift between mesh: " + impSlice + " and master mesh x,y,z: " +
		 (int)shift[0] + "," + (int)shift[1] + "," + (int)shift[2]);
		*/
		//modus..0: identify existing rings in master_mesh
		//modus..1: extend master_mesh with this meshes interior rings
		for(int modus = 0; modus < 2; ++modus) 
		{	
			int new_observations = 0;
			do
			{
				new_observations = 0;
				
				for(int i = 0; i < rings.size(); ++i)
				{
					Ring ringi = rings.get(i);
					if( (!ringi.is_interior) )
					{	continue;}
					if( (ringi.master != null) )
					{
						
						if( !validate_master_ring(ringi) )
						{
							System.out.println("Detected a corrupted ring - master pair");
							ringi.debug_info();
							ringi.master.debug_info();
							removefrom_master(); 
							//clean_up_master_mesh();
							return;
						}  
						continue;
					}
					int j0 = -1;
					for(int j = 0; (j < ringi.edges.size()) && (j0==-1); ++j)
					{
						final Bond bondj = ringi.edges.get(j);
						if(bondj.master != null)
						{	
							j0 = j;
							if( (bondj.a1.master == null) || (bondj.a2.master == null) )
							{
								System.out.println(" bond " + j0 + " has a master, but its ends dont");
								bondj.debug_info();
								bondj.master.debug_info();
								removefrom_master(); 
								return; 
							}
						} 
					}
					if(j0 == -1)
					{	continue;}
					Bond bondj = ringi.edges.get(j0); //must have master
					Atom a1 = bondj.a1; //must have master
					Atom a2 = bondj.a2; //must have master
					Ring ringj = bondj.get_other_ring(ringi); 
					
					Ring ringim = null;
					if(ringj.master != null)
					{
						ringim = bondj.master.get_other_ring(ringj.master); //may or may not exist
					}
					else
					{
						//ringi has an edge with a master, but the ring behind this edge does not
						//have a master, and therefore ringi should already have a master
						///TODO can we do something smarter about that
						System.out.println("unexpected topology");
						continue;	
					}
					
					
					//This triggers if something is going to be wrong however it cant be fixed
					if(ringim == null) //double check if we can do better than j0
					{
						for(int k = 0; k < ringi.edges.size(); ++k)
						{
							
							Bond bondk = ringi.edges.get(k); //must have master
							if(bondk.master != null)
							{
								Ring ringk = bondk.get_other_ring(ringi); //must have master
								if(ringk.master != null)
								{
									Ring ringkm = bondk.master.get_other_ring(ringk.master); //may or may not exist
									if(ringkm != null) //bondk works better than bondj, our earlier assumption was WRONG
									{
										j0 = k;
										bondj = bondk;
										a1 = bondk.a1;
										a2 = bondk.a2;
										ringj = ringk;
										ringim = ringkm;
										System.out.println("found better k than j0 when mastering ring: " + ringi.id + " in mesh: " +
											impSlice + " ringim: " + (ringim!=null?ringim.id:"null"));
									}
								}
							}
						
						}
					}

					if( (ringim != null) && (modus==0) ) //simple case just match up atoms and edges as needed
					{
						if( ringi.vertices.size() != ringim.vertices.size())
						{
							if(debug_mode)
							{
								System.out.println("conflicting topology in mesh: " + impSlice + " and master mesh");
								ringi.debug_info();
								ringim.debug_info();
								removefrom_master(); 
								throw new RuntimeException(Macro.MACRO_CANCELED);
							}
							ringi.is_interior = false;
							continue;
						}
						ringi.master = ringim;
						++ringim.observers;
						++new_observations;
						final int vs = ringi.vertices.size();
						
						Bond nb = a1.get_other_bond(ringi,bondj);
						Bond mnb = a1.master.get_other_bond(ringim,bondj.master);
						Atom na = nb.get_other_atom(a1);
						Atom mna = mnb.get_other_atom(a1.master);
						Atom ona = na;
						for(int v = 0; v < vs; ++v)//simply walk once around
						{
							if(mnb.left_ring == mnb.right_ring)
							{	System.out.println("Shit left and right ring are identical");} 
							if(nb.master == null)
							{
								nb.master = mnb;
								nb.master.observers++;	
							}
							else if(nb.master != mnb)
							{	System.out.println("Shit, nb.master != mnb");}
							
							if( (nb.master.left_ring != ringim) && (nb.master.right_ring != ringim) )
							{	System.out.println("Shit, nb.master does not touch ringim");}
							if(na.master == null)
							{
								na.master = mna;
								na.master.observers++;	
							}
							else if(na.master != mna)
							{	System.out.println("Shit, na.master != mna");}
							if(!na.master.rings.contains(ringim))
							{	System.out.println("Shit, na.master does not know ringim");}
							//push_atom_relations();
							nb = na.get_other_bond(ringi,nb);
							mnb = mna.get_other_bond(ringim,mnb);
							na = nb.get_other_atom(na);
							mna = mnb.get_other_atom(mna);	
						} 
						if(ona != na)
						{
							System.out.println("Not back at origin (" + (int)ona.pos[0] + "," + (int)ona.pos[1] +
							") after going n steps around n-ring!? actual end is (" + (int)na.pos[0] + "," + (int)na.pos[1]+")" );
						}
						//push_atom_relations();
						ringim.sort_all();
						
						if( !validate_master_ring(ringi) )
						{
							System.out.println("Ups, just assigned a non-fitting ring-master pair");
							ringi.debug_info();
							ringi.master.debug_info();
							bondj.debug_info();
							removefrom_master(); 
							//clean_up_master_mesh();
							return;
						}
									
					}
					else if( (ringim==null) && (modus==1) )//ringim does not yet exist in master_mesh
					{
						ringim = master_mesh.put_ring();
						ringi.master = ringim;
						ringim.observers = 1;
						++new_observations;
						ringim.is_closed = ringi.is_closed;
						ringim.is_interior = ringi.is_interior;
						ringim.pos = new float[3];
						ringim.pos[0] = ringi.pos[0]+shift[0];
						ringim.pos[1] = ringi.pos[1]+shift[1];
						ringim.pos[2] = ringi.pos[2]+shift[2];
						final int vs = ringi.vertices.size();
						
						Bond nb = bondj;
						Atom na = a1;
						
						Atom ona = na;
						for(int v = 0; v < vs; ++v,
													nb = na.get_other_bond(ringi,nb),
													na = nb.get_other_atom(na))//simply walk once around
						{
							Ring ringv = nb.get_other_ring(ringi);
							Ring ringvm = (ringv==null)?null:ringv.master;
							
							if(ringvm != null) //there is an mastered ring behind nb atoms and edges are known
							{
								if(ringvm == ringim)
								{	System.out.println("Shit ringvm is again ringim");}
								if( (nb.master == null) || (nb.a1.master == null) || (nb.a2.master == null) )
								{	System.out.println("A bond on a new master ring should have been known");}
								if(ringim.has_edge(nb.master) && 
									( (nb.master.left_ring == ringim)^(nb.master.right_ring == ringim) ) &&  
									(ringim.has_vertex(nb.master.a1) && ringim.has_vertex(nb.master.a2)) )
								{	
									System.out.println("A new edge was double checked");
									continue;
								}	
								
								if(nb.master.left_ring == nb.master.right_ring)
								{	System.out.println("An edge of a new master ring faces two times the same ring");}
								
								//All that we need to add to ringim should be already existing
								if(!ringim.neighbors.contains(ringvm))
								{
									ringim.neighbors.add(ringvm);
									ringvm.neighbors.add(ringim);
								}
								if(nb.master.right_ring != null)
								{	System.out.println("right_ring of an pre-existing bond.master is already occupied!?");}
								nb.master.right_ring = ringim;
							}
							else //the new bond does not directly connect somewhere in the master, ONE new atom is expected
							{	
								if(v==0)
								{
									System.out.println("Shit did not start start with bondj");
								}
								
								if(nb.master != null)
								{	System.out.println("Shit a new edge already exists bot did not connect to another master ring");}							
								
								if( ((nb.a1.master == null) && (nb.a2.master == null)))
								{	System.out.println("Shit the next bond along a new master ring has two loose ends");}
								
								if(na.master == null)
								{
									/* //This checks out
									if(na.rings.size() > na.edges.size())
									{
										System.out.println("Error: na has more rings than edges!");
										na.debug_info();
										//throw new RuntimeException(Macro.MACRO_CANCELED);
									}
									
									for(int r = 0; r < na.rings.size(); ++r)
									{	
										if( (na.rings.get(r).master != null) && (na.rings.get(r)!=ringi) )
										{	
											System.out.println("na.master is null but its ring: " +
											 na.rings.get(r).id + " has master: " + na.rings.get(r).master.id + "  v= " + v);
											
											validate_master_ring(na.rings.get(r));
											na.debug_info();
											na.rings.get(r).debug_info();
											na.rings.get(r).master.debug_info();
											//throw new RuntimeException(Macro.MACRO_CANCELED);
										}									
									}
									
									for(int r = 0; r < na.edges.size(); ++r)
									{
										if(na.edges.get(r).master != null)
										{	
											System.out.println("na.master is null but its bond: " +
											 na.edges.get(r).id + " has master: " + na.edges.get(r).master.id );
											na.debug_info();
											na.edges.get(r).debug_info();
											na.edges.get(r).master.debug_info();
											//throw new RuntimeException(Macro.MACRO_CANCELED); 
										}
									}
									*/
									
									float[] pos = new float[3];
									pos[0] = na.pos[0] + shift[0];
									pos[1] = na.pos[1] + shift[1];
									pos[2] = na.pos[2] + shift[2];
									
									na.master = master_mesh.put_atom(pos);
									na.master.observers = 1;
									na.master.is_interior = na.is_interior;
									na.master.shortest_distsqr = na.shortest_distsqr;
									
								}
								
								nb.master = master_mesh.put_bond();
								nb.master.observers = 1;
								nb.master.is_interior = nb.is_interior;
								nb.master.pos = new float[3];
								nb.master.pos[0] = nb.pos[0] + shift[0];
								nb.master.pos[1] = nb.pos[1] + shift[1];
								nb.master.pos[2] = nb.pos[2] + shift[2];
								nb.master.left_ring = ringim;
								nb.master.a1 = nb.a1.master;
								nb.master.a2 = nb.a2.master;
								nb.master.a1.edges.add(nb.master);
								nb.master.a2.edges.add(nb.master);
								nb.master.a1.neighbors.add(nb.master.a2);
								nb.master.a2.neighbors.add(nb.master.a1);
								nb.master.left_ring = ringim;	
							}
							
							ringim.edges.add(nb.master);
							ringim.vertices.add(na.master);
							na.master.rings.add(ringim);
							/* //This checks out
							if( (!(ringim.has_vertex(nb.master.a1))) && (!(ringim.has_vertex(nb.master.a2))) )
							{	System.out.println("Error nb.master.a1 and a2 are no vertices");}
							if(!(ringim.has_edge(nb.master)))
							{
								System.out.println("Error nb.master is no edge");
								if( (nb.master.left_ring == nb.master.right_ring) && (nb.master.left_ring != null) )
								{	System.out.println("Error* left and right ring of a nb.master are the same non nulls");}
								if( (nb.master.left_ring == ringim) || (nb.master.right_ring == ringim))
								{	System.out.println("Error** an yet not as edge registered bond does already know a ring");}
							}
							if( (nb.master.left_ring != ringim) && (nb.master.right_ring != ringim) )
							{	System.out.println("Error nb.master does not touch ringim");}
							*/	
						} 
						if(ona != na)
						{
							System.out.println("Not back at origin (" + (int)ona.pos[0] + "," + (int)ona.pos[1] +
							") after going n steps around a new n-ring!? actual end is (" + (int)na.pos[0] + "," + (int)na.pos[1]+")" );
						}
						//push_atom_relations();
						ringim.sort_all();
						
						if( !validate_master_ring(ringi) )
						{
							System.out.println("Ups just created a broken ring - master pair");
							ringi.debug_info();
							ringi.master.debug_info();
							removefrom_master(); 
							//clean_up_master_mesh();
							return;
						}												
					} //end if ringim && modus
					
				} //end for i in rings
			}
			while(new_observations > 0);
		} //end for modus 0 to 1
		
		
		for(int i = 0; i<rings.size(); ++ i)
		{
			Ring ringi = rings.get(i);
			if( (!ringi.is_interior) )
			{	continue;}
			if( (ringi.master != null) )
			{
				
				if( !validate_master_ring(ringi) )
				{
					System.out.println("Found a corrupted ring - master pair after adding all interior rings to master_mesh");
					ringi.debug_info();
					ringi.master.debug_info();
					removefrom_master(); 
					break;
				}  
				
				continue;
			}
			/* //Does happen if there are sepparate islands, not a real issue
			else
			{
				System.out.println("An interior ring was not added to the master_mesh");
				ringi.debug_info();
				removefrom_master();
				break; 
			}
			*/ 
		}
		System.out.println("added mesh: " + impSlice + " to master mesh");
		
	
	}
        
        public void manual_init()
        {
            
        }
	
	//no observations counting/mastering here, just push neighbor and edge relations to master_mesh
	private void push_atom_relations()  
	{
		for(int i = 0; i < atoms.size(); ++i)
		{
			final Atom atomi = atoms.get(i);
			if( (atomi.is_interior) && atomi.master!=null)
			{
				for(int k = 0; k < atomi.neighbors.size(); ++k)
				{
					Atom ak = atomi.neighbors.get(k);
					if( (ak.master != null) && (!atomi.master.neighbors.contains(ak.master)) )
					{
						atomi.master.neighbors.add(ak.master);
					}
				}	
				for(int k = 0; k < atomi.edges.size(); ++k)	
				{	
					Bond bk = atomi.edges.get(k);
					if( (bk.master != null) && (!atomi.master.edges.contains(bk.master)) )
					{
						atomi.master.edges.add(bk.master);
					}
				}	
				for(int k = 0; k < atomi.rings.size(); ++k)	
				{	
					Ring rk = atomi.rings.get(k);
					if( (rk.master != null) && (!atomi.master.rings.contains(rk.master)) )
					{
						atomi.master.rings.add(rk.master);
					}		
				}
			}
		}
		
	}
	
	
	
	private void removefrom_master()
	{
		//this has also be done for detached atoms
		for(int i=0; i < atoms.size(); ++i)
		{
			Atom ma = atoms.get(i).master;
			if( (ma != null) && master_mesh.atoms.contains(ma) && (!ma.detached) )
			{
				atoms.get(i).master = null;
				if( (ma.observers > 0) && (--ma.observers < 1) )
				{	ma.detach();}
			}
		}
		//In thoery detaching all atoms should fully suffice
		//TODO test if the other two loops can be omitted
		for(int i=0; i < bonds.size(); ++i)
		{
			Bond ma = bonds.get(i).master;
			if( (ma != null) && master_mesh.bonds.contains(ma) && (!ma.detached) )
			{
				bonds.get(i).master = null;
				if( (ma.observers > 0) && (--ma.observers < 1) )
				{	ma.detach();}
			}
		}	
		for(int i=0; i < rings.size(); ++i)
		{
			Ring ma = rings.get(i).master;
			if( (ma != null) && master_mesh.atoms.contains(ma) && (!ma.detached) )
			{
				rings.get(i).master = null;
				if( (ma.observers > 0) && (--ma.observers < 1) )
				{	ma.detach();}
			}
		}		
		clean_up_master_mesh();
	}
	
	private boolean validate_master_ring(Ring ringi)
	{
		if(ringi.master == null)
		{
			System.out.println("ringi.master is null");
			return false;
		}
		
		if( (ringi.vertices.size() != ringi.master.vertices.size()) ||
			(ringi.edges.size() != ringi.master.edges.size()) )
		{
			System.out.println("number of edges/vertices dont match");
			return false;
		}
		for(int j = 0; j < ringi.vertices.size(); ++j)
		{
			Atom atomj = ringi.vertices.get(j);
			Bond bondj = ringi.edges.get(j);
			
			if( (atomj.master == null) || (bondj.master == null) )
			{	
				System.out.println("vertex or edge "+ j + " dont have a master");
				if((atomj.master == null))
				{	System.out.println("vertex.master is null");}
				if((bondj.master == null))
				{	System.out.println("edge.master is null");}
				return false;
			}
			
			if( !( (ringi.master == bondj.master.left_ring)^
				   (ringi.master == bondj.master.right_ring) ) )
			{	
				System.out.println("ringi.master is not knowm by bondj.master, j: " + j);
				bondj.debug_info();
				bondj.master.debug_info();
				return false;
			}
			
			
			
			if( (!ringi.master.has_vertex(atomj.master)) ||
				(!ringi.master.has_edge(bondj.master)) )
			{	
				System.out.println("ringi.master does not know atomj.master or bondj.master: " + j);
				return false;
			}
			
			if (!atomj.master.rings.contains(ringi.master) )
			{	
				System.out.println("atomj.master does not know ringi.master j: " +
				 j + " x,y: " + (int)atomj.pos[0] + "," + (int)atomj.pos[1] );
				return false;
			}		
			
			for(int k = 0; k < atomj.edges.size(); ++k)
			{
				Bond bondk = atomj.edges.get(k);
				if(ringi.has_edge(bondk))
				{
					if ( !atomj.master.edges.contains(bondk.master) ) 
					{
						System.out.println("bondk.master is " + ((bondk.master==null)?"null":"unkown") );	
						if(bondk.master != null)
						{
							System.out.println("bondk id: " + bondk.id + " x,y: " +
								(int)bondk.pos[0] + "," + (int)bondk.pos[1]);
							System.out.println("bondk.master id: " + bondk.master.id + " x,y: " +
								(int)bondk.master.pos[0] + "," + (int)bondk.master.pos[1]);
						}
						
						System.out.println("atomj.master.edges.size(): " + atomj.master.edges.size());
						for(int i = 0; i < atomj.master.edges.size(); ++ i)
						{
							System.out.println("atomj.master.edges.get("+i+").id: " + atomj.master.edges.get(i).id + 
							" x,y: " + (int)atomj.master.edges.get(i).pos[0] + "," + (int)atomj.master.edges.get(i).pos[1] );
						}
						System.out.println("atomj.edges.size(): " + atomj.edges.size());
						for(int i = 0; i < atomj.edges.size(); ++ i)
						{
							System.out.println("atomj.edges.get("+i+").id: " + atomj.edges.get(i).id + 
							" x,y: " + (int)atomj.edges.get(i).pos[0] + "," + (int)atomj.edges.get(i).pos[1] );
							if(atomj.edges.get(i).master != null)
							{
								System.out.println("atomj.edges.get("+i+").master.id: " + atomj.edges.get(i).master.id + 
								" x,y: " + (int)atomj.edges.get(i).master.pos[0] + "," + (int)atomj.edges.get(i).master.pos[1] );
							}
							else
							{
								System.out.println("atomj.edges.get("+i+").master is null");
							}
						}
					
					
					}
				}
			}
			
			
			if( !( (bondj.a1.master == bondj.master.a1)^(bondj.a1.master == bondj.master.a2 ) &&
				   (bondj.a2.master == bondj.master.a1)^(bondj.a2.master == bondj.master.a2 )	  )
			 )
			{	
				System.out.println("atoms and their masters dont pair up for bondj and bondj.master, j: " + j);
				bondj.debug_info();
				bondj.master.debug_info();
				System.out.print('\n');
				
				bondj.a1.debug_info();
				System.out.print('\n');
				bondj.a1.master.debug_info();
				System.out.print('\n');
				bondj.master.a1.debug_info();
				
				System.out.print('\n');
				bondj.a2.debug_info();
				System.out.print('\n');
				bondj.a2.master.debug_info();
				System.out.print('\n');
				bondj.master.a2.debug_info();
				
				System.out.print('\n');
				return false;
			}
			
			
			if(bondj.master.left_ring == bondj.master.right_ring)
			{	
				System.out.println("bondj master is touching ring master twice");
				bondj.debug_info();
				bondj.master.debug_info();
				return false;
			} 
			
			
		}	
		return true;				
	}

}
