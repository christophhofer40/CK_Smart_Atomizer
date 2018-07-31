import ij.*;
import ij.io.*;
import ij.io.FileInfo;
import ij.process.*;
import ij.gui.*;

import ij.plugin.filter.*;
import ij.gui.DialogListener;
import ij.plugin.filter.PlugInFilterRunner;

import java.awt.*;
import java.awt.event.*;
import java.awt.Checkbox;
import java.util.Arrays;

import java.io.*;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.FileReader;

import krami.SurfaceMesh;
import krami.SurfaceMesh.*;
import krami.HexPixels;

class AtomJoinerCore implements DialogListener, MouseListener
{
	int flags;
	int DONE = 4096;
	public ImagePlus imp = null;
	ImageStack impSt = null;
	int pwidth = 2;
	int pheight = 1;
	final int panelsep = 5;
	final float radius = 5;
	ImagePlus panels = null;
	ImageProcessor panelP = null; 
	int psl1 = 1;
	int psl2 = 2;
	int sl1 = psl1;
	int sl2 = psl2;
	
	File topFile,outFile;
	static SurfaceMesh[] meshes = null;
	static SurfaceMesh master_mesh = null;
	NonBlockingGenericDialog gd = null;
	boolean hex_pixels = false;
	HexPixels imgHP = null;
	boolean periodic = false;
	boolean inverted = false;
	boolean reset = false;
	boolean reset_position = false;
	boolean place_atoms = false;
	int impWidth = -1;
	int impHeight = -1;
	
	Roi lroi = null;
	
    // these listeners are activated if the selection is changed in the corresponding ImagePlus
    public void mousePressed(MouseEvent e) {};//{IJ.log("pressed");}   
    public void mouseDragged(MouseEvent e) {};//{IJ.log("drag");}
    public void mouseClicked(MouseEvent e) 
    { 
		Roi nroi = panels.getRoi();
		if(nroi==null || nroi.getType()!=Roi.POINT)
		{	run(null);} 
	};//redraw atoms and cancel lines
    public void mouseReleased(MouseEvent e) //{IJ.log("release");}
    {
		Roi nroi = panels.getRoi();
		if( (nroi != null) && (gd != null) && (gd.isPreviewActive()))
		{
			if(lroi == null || nroi != lroi)
			{
				lroi = nroi;
				switch(nroi.getType())
				{
					case Roi.LINE:
						process_Line(nroi);
					break;
					case Roi.OVAL:
						process_Circle(nroi);
					break;
					case Roi.RECTANGLE:
						process_Box(nroi);
					break;
					case Roi.POINT:
						process_Point(nroi);
					default:
				}
				
				
				
			}
		}
		
	};
    public void mouseExited(MouseEvent e) {};//{	IJ.log("mouse left");}
    public void mouseEntered(MouseEvent e) {};//{ IJ.log("mouse entered");}
    public void mouseMoved(MouseEvent e) {}
   
	
	
	
	
	AtomJoinerCore(int flg) //dont call that with flg == DONE == 4096
	{
		flags = flg;
	};

	public int showDialog(final ImagePlus imp, final String command, final PlugInFilterRunner pfr)
	{
		gd = new NonBlockingGenericDialog( command ); 
		gd.addMessage( imp.getTitle() );
		gd.addMessage( panels.getTitle() );
		gd.addNumericField("left panel: ", psl1, 0 );
		gd.addNumericField("right panel: ", psl2, 0 );
		gd.addCheckbox("put/nudge atom", place_atoms);
		gd.addCheckbox("reload", reset);
		if(pfr != null)
		{	
			gd.addPreviewCheckbox(pfr,"interactive");
			gd.addDialogListener(this);	
		}
		else
		{	gd.addMessage("preview unavailable");}
		
		gd.showDialog();
		if(panels!=null)
		{
			ImageWindow win = panels.getWindow();
			ImageCanvas canvas = win.getCanvas();
			canvas.removeMouseListener(this);
        }
        if (gd.wasCanceled())
		{	return DONE;}
		writetopfile();
		return flags;
	}
	
	
	public boolean dialogItemChanged(final GenericDialog gd, final AWTEvent e) 
	{
		boolean valid_input = true;
		sl1 = (int)gd.getNextNumber();
		sl2 = (int)gd.getNextNumber();
		valid_input &= (sl1 > 0);
		valid_input &= (sl2 > 0);
		valid_input &= (sl1 != sl2);
		place_atoms = gd.getNextBoolean();
		reset = gd.getNextBoolean();
		//valid_input &= (sl1 != psl1);
		//valid_input &= (sl2 != psl2);
		valid_input &= !gd.invalidNumber();
		return ( valid_input );
	
	}
	
	public int init(ImagePlus nimp)
	{
		boolean keep_going = true;
		imp = nimp;
		FileInfo fi = imp.getOriginalFileInfo();
		String path = (fi!=null)?fi.directory:IJ.getDirectory("current");
		String topname = imp.getShortTitle() + ".top";
		OpenDialog od = new OpenDialog("Open top file",path,topname);
		String topT = od.getPath();
		if(topT != null)
		{
			topFile = new File(topT);
			if( (topFile.isFile() && topFile.canRead()))
			{	IJ.log(topT);}
			else
			{	return 1;}
		}
		else
		{	
			keep_going = IJ.showMessageWithCancel("Manual Mode","Proceed without any inital atoms?");
			reset = false;
		}
		if(!keep_going)
		{	return 1;}
		if(reset)
		{
			reset_position = IJ.showMessageWithCancel("reseting coordinates","reset master positions to first observations");  	
		}
		
		
		impWidth = imp.getWidth();
		impHeight = imp.getHeight();
		
		pwidth = 2*impWidth + panelsep;
		pheight = impHeight;
		//IJ.log("pwidth: " + pwidth + "  height: " + pheight );
		impSt = imp.getStack();
		int bitDepth = impSt.getBitDepth();
		int impDepth = impSt.getSize();
		if(panels == null)
		{
			ImageStack panelSt = ImageStack.create(pwidth,pheight,1,bitDepth);
			panels = new ImagePlus("panels", panelSt);
		}
		panelP = panels.getProcessor();
		
		update_panel();
		
		
		
		
		
		SurfaceMesh.hex_pixels_default = hex_pixels;
		SurfaceMesh.periodic_default = periodic; 
		SurfaceMesh.inverted_default = inverted; 
			
		meshes = new SurfaceMesh[impDepth];	
			
		for(int i = 0; i < impDepth; ++i)
		{
			meshes[i] = new SurfaceMesh(impWidth,impHeight,i+1);
			meshes[i].hex_pixels = hex_pixels;
			meshes[i].imgHP = imgHP;
			meshes[i].periodic = periodic;
			meshes[i].inverted = inverted;
		}
		AtomJoinerCore.master_mesh = new SurfaceMesh(impWidth,impHeight,0);
		SurfaceMesh.master_mesh = AtomJoinerCore.master_mesh;
		master_mesh.hex_pixels = hex_pixels;
		master_mesh.imgHP = imgHP;
		master_mesh.periodic = periodic;
		master_mesh.inverted = inverted;
		if(topT != null)
		{	reatopfile();}
		
		panels.show();
		
		ImageWindow win = panels.getWindow();
        ImageCanvas canvas = win.getCanvas();
        canvas.addMouseListener(this);
	
		return 0;
	}
	
	void update_panel()
	{
		ImageProcessor aSlice = impSt.getProcessor(sl1);
		panelP.insert(aSlice, 0, 0);
		ImageProcessor bSlice = impSt.getProcessor(sl2);
		panelP.insert(bSlice, impWidth + panelsep, 0);
		psl1 = sl1;
		psl2 = sl2;
		if(panels.isVisible());
		{	panels.updateAndRepaintWindow();}
		//IJ.log("new panels, left: " + psl1 + "  right:" + psl2);
	}
	
	
	public void process_Line(Roi roi)
	{
		//IJ.log("processing Line");
		Polygon pts = roi.getPolygon();
		int x1 = (int) pts.xpoints[0];
		int y1 = (int) pts.ypoints[0];
		int x2 = (int) pts.xpoints[1];
		int y2 = (int) pts.ypoints[1];
		if(x1 > impWidth)
		{
			int xt = x1;
			x1 = x2; 
			x2 = xt;
			int yt = y1;
			y1 = y2; 
			y2 = yt;
		}
		if(x1 < impWidth && x2 >= impWidth + panelsep)
		{
			x2-= (impWidth+panelsep);
			//IJ.log("(" + x1 + "," + y1 + ") -> (" + x2 + "," + y2 + ")");
			float[] p1 = {(float)x1,(float)y1,0.0f};
			float[] p2 = {(float)x2,(float)y2,0.0f};
			
			Atom[] at1 = meshes[psl1-1].findAtoms(p1,radius);
			Atom[] at2 = meshes[psl2-1].findAtoms(p2,radius);
			
			if( (at1 != null) && (at1.length > 0) && (at2 != null) && (at2.length > 0) )
			{
				Atom a1 = at1[0];
				Atom a2 = at2[0];
				
				if( (a1.master != null) && (a2.master != null) )
				{
					if(a1.master == a2.master)
					{
						//IJ.log("canceling an identfication");
						if( (a1.master.observers-=2) < 1)
						{	a1.master.detach();}
						a1.master = null;
						a2.master = null;
						a1.is_interior = false;
						a2.is_interior = false;
					}
					else
					{
						//IJ.log("These atoms are different");
					}	
				}
				else if ((a1.master == null) && (a2.master == null))
				{
					//IJ.log("creating fresh master atom");
					Atom ma = master_mesh.put_atom(Arrays.copyOf(a1.pos,a1.pos.length));
					a1.master = ma;
					a2.master = ma;
					ma.observers+=2;
					a1.is_interior=true;
					a2.is_interior=true;
				}
				else
				{
					//IJ.log("adding observers to existing master");
					Atom ma = (a1.master==null)?a2.master:a1.master;
					final int last_observers = ma.observers;
					if(a1.master==null) 
					{
						Atom[] observers = meshes[psl1-1].get_observers(ma);
						if( (observers==null) || (observers.length == 0) )
						{
							a1.master = ma;
							++(ma.observers);
							a1.is_interior = true;	
						}
						else
						IJ.log("there could only be one unique observation in the left view");
						
					}
					//else	
					if(a2.master==null)
					{
						Atom[] observers = meshes[psl2-1].get_observers(ma);
						if( (observers==null) || (observers.length == 0) )
						{
							a2.master = ma;
							++(ma.observers);
							a2.is_interior = true;
						}
						else
						IJ.log("there could only be one unique observation in the right view");
					}
					//if(ma.observers != last_observers + 1)
					//{	IJ.log("UPS messed up linking another atom to master mesh");}
				
				}
			}
			else
			{
				//IJ.log("no nearby Atoms");
			}
			
			
		}
		run(null);		
	}
	
	
	public void process_Box(Roi roi)
	{
		
		Rectangle  rect = roi.getBounds();
		int x1 = (int) rect.x;
		int y1 = (int) rect.y;
		int x2 = (int) x1+rect.width;
		int y2 = (int) y1+rect.height;
		if(! ((x2 < impWidth) || (x1 >= impWidth + panelsep) ))
		{	return;}
		//IJ.log("Processing Box");
		int src = psl1-1;
		int dst = psl2-1;
		boolean right_side = false;
		final int pshift = impWidth + panelsep;
		if(x1 >= pshift)
		{
			right_side = true;
			src = psl2-1;
			dst = psl1-1;
			x1 -= pshift;
			x2 -= pshift;
		}
		
		
		
		Atom[] boxed = meshes[src].boxed_atoms((float)x1,(float)y1,(float)x2,(float)y2);
		//if(boxed == null || boxed.length==0)
		//{	IJ.log("Empty box");}
		//IJ.log("src: " + (src+1) + " dst: " + (dst+1));
		Overlay ov = panels.getOverlay();
		for(int i = 0; i < boxed.length; ++i)
		{
			Atom atomi = boxed[i];
			Atom ma = atomi.master;
			if(ma == null) continue;
			Atom[] others = meshes[dst].get_observers(ma);
			if( (others == null) || (others.length == 0) ) continue; 
			Atom atomj = others[0];
			Line line = new Line((int)atomi.pos[0] + (right_side?pshift:0),(int)atomi.pos[1],
								(int)atomj.pos[0] + (right_side?0:pshift),(int)atomj.pos[1]);
			line.setStrokeColor(Toolbar.getForegroundColor());
			ov.add(line);					
		}
		panels.setOverlay(ov);
			
	}
	
	public void process_Circle(Roi roi)
	{
		//IJ.log("Processing Circle");
		Rectangle  rect = roi.getBounds();
		float x1 = (float) (rect.x+0.5f*rect.width);
		float y1 = (float) (rect.y + 0.5f*rect.height);
		float z1 = 0.0f;
		boolean right_side = true;
		int sl = -1;
		if(x1 < impWidth)
		{	sl = psl1;}
		else if (x1 >= impWidth+panelsep)
		{	
			sl = psl2;
			right_side = false;
			x1-=(impWidth+panelsep);
		}
		
		float[] p1 = {x1,y1,z1};
		float r = (float)Math.sqrt(0.25*(rect.width*rect.width+rect.height*rect.height));
		Atom[] pa = meshes[sl-1].findAtoms(p1,r);
		Atom atom = null;
		if(pa == null || pa.length == 0)//simply put an atom
		{	
			if(place_atoms)
			{	meshes[sl-1].put_atom(p1);}
		}
		else
		{
			if(place_atoms) //hijack and spare the closest atoms
			{	pa[0].pos = p1;} 
			for(int i = (place_atoms?1:0); i < pa.length; ++i)
			{
				Atom atomi = pa[i];
				if(atomi.master != null)
				{
					Atom am = atomi.master;
					--am.observers;
					if(am.observers == 0)
					{	
						am.is_interior = false;
						am.detach();
					}
					atomi.master = null;
				}
				atomi.is_interior = false;
				atomi.detach();
			}
		
		}
		run(null);
			
	}
	
	public void process_Point(Roi roi)
	{
		Rectangle  rect = roi.getBounds();
		float x1 = (float) (rect.x);
		float y1 = (float) (rect.y);
		float z1 = 0.0f;
		boolean left_side = true;
		int sl = -1;
		if(x1 < impWidth)
		{	sl = psl1;}
		else if (x1 >= impWidth+panelsep)
		{	
			sl = psl2;
			left_side = false;
			x1-=(impWidth+panelsep);
		}
		
		float[] p1 = {x1,y1,z1};
		Atom[] pa = meshes[sl-1].findAtoms(p1,2*radius);
		if(pa != null && pa.length > 0)
		{
			Atom atomi = pa[0];
			final float x2 = left_side?atomi.pos[0]:atomi.pos[0]+impWidth+panelsep; 
			final float y2 = atomi.pos[1]; 
			atomi.debug_info();
			if(atomi.master != null)
			{	atomi.master.debug_info();}
			Overlay ov = panels.getOverlay();
			
			OvalRoi nr = new OvalRoi((int)(x2-2*radius),(int)(y2-2*radius),(int)(4*radius),(int)(4*radius));
			nr.setStrokeColor(Toolbar.getForegroundColor());
			ov.add(nr);
			panels.setOverlay(ov);	
		}
	
	
	}
	
	public void run(ImageProcessor ip) 
	{
		if(reset)
		{	
			init(imp);
			reset = false;
			if(gd!=null)
			{	
				Checkbox cb = (Checkbox)gd.getCheckboxes().get(1);
				cb.setState(false);	
			}
		}
		
		
		
		ImageProcessor ipp = panels.getProcessor();
		//IJ.log("sl1: " + sl1 + "  sl2:" + sl2);
		if( (psl1 != sl1) || (psl2 != sl2) )
		{	update_panel();}
		
		
		master_mesh.refresh();
			
		meshes[sl1-1].fake_refresh();
		meshes[sl1-1].drawMarks(imp, true, false, false, false);
		
		meshes[sl2-1].fake_refresh();
		meshes[sl2-1].drawMarks(imp, true, false, false, false);
		
		Overlay ov = imp.getOverlay();
		Roi[] rois = ov.toArray();
		Overlay po = new Overlay();
		
		if( (rois != null) && (rois.length>0) )
		{
			for(int i = 0; i < rois.length; ++i)
			{
				Roi roi = rois[i];
				final int sl = roi.getPosition();
				if( (sl!=sl1) && (sl!=sl2) ) continue;
				Roi cp = (Roi)roi.clone();
				cp.setImage(panels);	
				cp.setPosition(1);
				if(sl==sl2)
				{	
					int x = cp.getBounds().x + impWidth + panelsep;
					int y = cp.getBounds().y;
					cp.setLocation(x,y);
				}
				po.add(cp);	
			}	
			
		}
		Overlay oo = panels.getOverlay();
		if(oo != null ) oo.clear();
		panels.setOverlay(po);
		return;	
	}
	
	private void reatopfile()
	{
		IJ.log("reading " + topFile.getPath());
		try
		{
			int slice = -1;//invalid number
			BufferedReader reader = new BufferedReader(new FileReader(topFile));
			
			while (reader.ready())
			{
				String line = reader.readLine();
				if(line == null)
				{	break;}
				if(line.startsWith("#")) continue; //skip comments
				String[] words = line.split("\\s+");
				if(words==null) continue; //skip lines of only blanks
				if(words[0].equals("MASTER") && slice == -1)
				{
					slice = 0; //also invalid for invisible master
				}
				else if(words[0].equals("VIEW") )
				{
					slice = Integer.parseInt(words[1])+1;					 
				}
				else if(words[0].equals("ATOM"))
				{
					SurfaceMesh mesh = (slice==0)?master_mesh:meshes[slice-1];
					int id = Integer.parseInt(words[1]);
					float[] pos = new float[]{0.0f,0.0f,0.0f};
					int s = (slice==0?1:0);
					pos[0] = Float.parseFloat(words[2+s]);
					pos[1] = Float.parseFloat(words[3+s]);
					//pos[2] = Float.parseFloat(words[4+s]); //maybe problematic with mouse clicking 
					
					Atom ma = mesh.put_atom( pos );
					ma.intensity = 1.0f;
					if( (id != -1) && (mesh!=master_mesh) )
					{
						//the master may be actually missing
						ma.master = master_mesh.getAtom(id);
						if(ma.master!= null)
						{
							++(ma.master.observers);
							ma.is_interior = true;
							if(reset_position && ma.master.observers==1)
							{
								ma.master.pos[0] = ma.pos[0];
								ma.master.pos[1] = ma.pos[1];
								ma.master.pos[2] = 0.0f;
							}
						}
					} 	
				}
				//ignore BONDS,RINGS,SKEWMATRIX,QUATERNION, and whatever else
				
			}
			reset_position = false;
			master_mesh.refresh();
			for(int i=0; i < imp.getStackSize(); ++i )
			{
				meshes[i].fake_refresh();
				meshes[i].drawMarks(imp, true, false, false, false);
			}
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	
	private void writetopfile()
	{
		//outFile = new File(topFile.getPath()+".out");
		master_mesh.refresh();
		for(int v = 0; v<meshes.length; ++v)
		{
			meshes[v].keep_atom_masters=true;
			meshes[v].refresh();
			meshes[v].drawMarks(imp, true, true, true, false);
			meshes[v].addto_master(true);
			meshes[v].keep_atom_masters=false;
		}
		
		//master_mesh.refresh();
		outFile = new File(topFile.getPath());
		IJ.log("rewriting " + outFile.getPath());
	
		try
		{
			FileWriter fw = new FileWriter(outFile, false);
			PrintWriter pw = new PrintWriter(fw);
			
			pw.println("#creator: " + getClass().getSimpleName() + " source: " + imp.getTitle());
			pw.println("#MASTER are the ATOMs in the joined topology with their id, view counts and approximate x,y,z coordinates\n" +
					   "#BONDs connect the two ATOMs with id1 and id2\n" +
					   "#RINGs list the ordered ATOM ids of their vertices" +	
					   "#Every VIEW lists a LABEL, a TRANSLATION, a SKEWMATRIX,\n" +
					   "#a rotation QUATERNION as well as all ATOMS with their id, x and y in image pixels and z=0.0");
			pw.println("\n\nMASTER\t" + master_mesh.fatoms.length);
			for(int i = 0; i < master_mesh.fatoms.length; ++i)
			{
				final Atom am = master_mesh.fatoms[i];
				pw.println("ATOM\t" + am.id + "\t" + am.observers + "\t" + am.pos[0] + "\t" + am.pos[1] + "\t" + am.pos[2]);
			}		   
			pw.print('\n');
			for(int i = 0; i < master_mesh.fbonds.length; ++i)
			{
				final Bond bondi = master_mesh.fbonds[i];
				/*
				if( (bondi.left_ring == null) || (bondi.right_ring == null) ||
				    (!bondi.left_ring.is_interior) || (!bondi.right_ring.is_interior) )
				{	continue;}
				*/
				pw.println("BOND\t" + bondi.a1.id + "\t" + bondi.a2.id);
			}
			pw.print('\n');
			for(int i = 0; i < master_mesh.frings.length; ++i)
			{
				final Ring ringi = master_mesh.frings[i];
				pw.print("RING");
				for(int j = 0; j < ringi.vertices.size(); ++j)
				{	pw.print("\t"+ ringi.vertices.get(j).id);}
				pw.print('\n');
			}				
			
			
			for(int v = 0; v<meshes.length; ++v)
			{
				if( (meshes[v].fatoms == null) || (meshes[v].fatoms.length == 0) )
				{	continue;}
				String lbl = impSt.getSliceLabel(v+1);
				lbl = (lbl!=null)?lbl.split("\\n")[0]:"null";
				pw.println("\n\nVIEW\t" + v);
				pw.println("#view of slice:" + (v+1) + " in " + imp.getTitle());
				pw.println("LABEL\t"+lbl);
				pw.println("TRANSLATION\t0.0\t0.0\t0.0");
				pw.println("SKEWMATRIX\t1.0\t0.0\t0.0\t1.0");
				pw.println("QUATERNION\t0.0\t0.0\n");
				
				for(int i = 0; i < meshes[v].fatoms.length; ++i)
				{
					final Atom a = meshes[v].fatoms[i];
					final Atom am = a.master;
					if(am==null)
					{	
						pw.println("ATOM\t" + (-1) + "\t" + a.pos[0] + "\t" + a.pos[1] + "\t" + a.pos[2]);	
					}
					else
					{
						pw.println("ATOM\t" + am.id + "\t" + a.pos[0] + "\t" + a.pos[1] + "\t" + a.pos[2]);
					}
				}	
			}
			pw.close();
			IJ.log("rewrote " + outFile.getName());
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	
	}

}
