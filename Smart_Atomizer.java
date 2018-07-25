import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.*;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.plugin.*;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.gui.Overlay;
import ij.plugin.frame.RoiManager;
import java.awt.*;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Choice;
import java.awt.Checkbox;
import java.awt.TextField;
import java.util.Random;
import java.util.Arrays;
import java.util.Scanner;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.LinkedList;
import java.io.*; //File BufferedReader Reader etc.


import krami.SurfaceMesh;
import krami.SurfaceMesh.*;
import krami.HexPixels;
import krami.SharedMerit;


import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.FileReader;

public class Smart_Atomizer implements PlugIn
{
	public static ImageStack next_modelSt = null;
	public static boolean external_merit = false;
	public static String external_merit_source = null;
	public static SharedMerit sharedmerit = null; //ought to be set by Hex_Magic
	private static LinkedList<NonBlockingGenericDialog> dialogs = new LinkedList<NonBlockingGenericDialog>();
	//these may be accessed and modified by any thread	
	
	static SurfaceMesh[] meshes = null;
	SurfaceMesh mesh = null;
	static SurfaceMesh master_mesh = null;
	static String master_sig_string = null;
	
	static HexPixels atomHP = null;
	static HexPixels imgHP = null;
	static HexPixels tinyHP = null;
	static HexPixels boxHP = null;
	
	static public int get_instance_count()
	{	return instance_count.get();}
	
	static AtomicInteger instance_count = new AtomicInteger(0);
	static AtomicInteger nextHelperSlice = new AtomicInteger(0);
	static AtomicInteger completed_cycles = new AtomicInteger(0);
	static AtomicIntegerArray taken_slices = new AtomicIntegerArray(16);
	
	//static volatile is shared and may be read and set by anyone, 
	//but no self-conditional changes like &=,|=, ...  
	static volatile boolean reset_enhanced = false;
	
	//static "non-volatile" members are only allowed to be modified by the master
	//when there are no slaves mostly do_Dialog() and validate_Input() and special cases
	//in run() and run_slice(), every thread can see them
	
	static boolean inverted = SurfaceMesh.inverted_default; //false;
	static boolean periodic = SurfaceMesh.periodic_default; //false;
	static boolean hex_pixels = SurfaceMesh.hex_pixels_default; //false
	static boolean do_stack = false;
	static boolean do_setup_atom = false;
	static boolean enable_holes = false;
	static boolean enforce_solid = false;
	static boolean equalize_intensities = false;
	static boolean optimize_enhanced = false;
	static boolean protect_atoms = false;
	static boolean fix_topology = false;
	static boolean show_topology = false;
	static boolean optimize_fine = false;
	static boolean optimize_weights = true;
	static boolean optimize_FWHM = false;
	static boolean optimize_Brightness = true;
	static boolean count_atoms = false;
	static boolean read_sig_file = false;
	static boolean join_meshes = false;
	static boolean write_sig_file = false;
	static boolean write_top_file = false;
	static boolean append_another_task = false;
	static boolean mesh_from_xyz = false;
	
	static Random random = new Random(System.currentTimeMillis());
	static ImagePlus imp = null;
	static ImagePlus atom = null;
	static ImagePlus simulated = null;
	static ImagePlus enhanced = null;
	static ImagePlus dual = null;
	static ImagePlus deviations = null;
	
	static String impT = "";
	static String atomT = "atom.tif";
	static String simulatedT = "sim.tif";
	static String enhancedT = "model.tif";
	static String deviationsT = "diff.tif";
	
	static int impWidth = -1;
	static int impHeight = -1;
	static int impDepth = 0;
	static int impArea = -1;
	private int validArea = -1;
	
	
	static int impSlice0 = 1;
	static int cycles = 1;
	static int grand_runs = 2000;
	static double wobble = 2.0;
	static int fine_runs = 2;
	static int atomWidth = -1;
	static int atomHeight = -1;
	static double mix = 1.0;
	static double atomFWHM = 6;
	static double atomSolidD0 = 6;
	static double peakWeightMin = 0.1;
	
	static int numTasks = 0;
	static int maxThreads = 8;
	static Point loc = null;
	static Rectangle rec = null;
	
	static float[] a_pix = null;
	
	//ImageStack appear to be threadsafe, every thread shall
	//only access his reserved frame with pull_pixels() and push_pixels()
	
	static ImageStack impSt = null;
	static ImageStack atomSt = null;
	static ImageStack simulatedSt = null;
	static ImageStack enhancedSt = null;
	static ImageStack dualSt = null;
	static ImageStack deviationsSt = null;
	
	//regular globals are per Thread and need to be initalized
	//mostly ensural_all_images() and init_blabla() therin 
	
	boolean waited = false;
	
	float[] imp_pix = null;
	float[] sim_pix = null;
	float[] enh_pix = null;
	float[] dual_pix = null;
	float[] dev_pix = null;
	
	boolean very_first_slice = false;
	int impSlice = 1;
	int seat = 0;
	
	double impTotal = 0.0;
	double impAvg = 0.0;
	double impStd = 0.0;
	double model_merit = 0.0;
	
	double total_peak_weight = 0.0;
	double total_hole_weight = 0.0;
	
	double avg_peak_weight = 0.0;
	double avg_hole_weight = 0.0;
	int numAtoms = 0;
	int numHoles = 0;
	int deltaAtoms = 0;
	
	double atomSolidD = 7.5;
	double modelFloor = 0.0;
	
	TextRoi status_roi = null;
	
	public void run(String arg) 
	{
		try
		{
			seat = instance_count.getAndIncrement();
			if (seat > 0)
			{	
				if( (taken_slices != null) && (seat < taken_slices.length()) && (seat < numTasks) )
				{
					//System.out.println(getClass().getSimpleName()+": slave is joining, thread group size: " +
					//					instance_count.get() + " task slots: " + numTasks);
				}
				else
				{
					seat = instance_count.decrementAndGet();
					System.out.println( getClass().getSimpleName() + ": thread: "+ (seat+1) +
								" has no task, thread slots: " + numTasks +
								 " nextSlice: " + nextHelperSlice.get() +
								  " cyles: " + completed_cycles.get() + "/" + cycles +
								  " remaining tasks: " + dialogs.size() +
								   " running: " + instance_count.get() +
								   " busy slices: " + taken_slices	);
					if(IJ.showMessageWithCancel(getClass().getSimpleName(), "Dismiss all further runs & tasks?"))
					{
						System.out.println("Dismissing all further runs and tasks for " + getClass().getSimpleName() );
						dialogs.clear();
						nextHelperSlice.set(0);
					}
					return;
				}
				
			}
			else if(seat == 0)
			{
				//System.out.println(getClass().getSimpleName()+": master is starting, thread pool size: " + instance_count.get() );
				if( (imp == null) || (!imp.isVisible()) )
				{
					imp = IJ.getImage();
					impT = imp.getTitle();
				}
				/*if(imp == null)
				{
					IJ.noImage();
					seat = instance_count.decrementAndGet();
					return;
				}*/
				dialogs.clear();
			}
			else //really bad error seat < 0
			{
				System.out.println(getClass().getSimpleName()+": threading ERROR, seat: " + seat + ", instances: " + instance_count.get());
				System.out.println("resetting both to zero, you can use Thread_Killer.java to kill zombie threads");
				instance_count.set(0);
				return;
			}
			
			do
			{
				if(seat == 0)
				{
					numTasks = 1; //make sure no slave can join before we need him
					instance_count.set(1); //undo the false decrement for the very last thread
					nextHelperSlice.set(0);// any slave that has sneak in will be told to quit
				}
				if( (seat != 0) || doDialog() )
				{
					if( seat == 0 ) //we are the master
					{
						if( validateInput(true) != 0 ) //true = init_stuff
						{	continue;} //the master will quit if there are still other threads
						completed_cycles.set(0);
						nextHelperSlice.set(impSlice);
						if(sharedmerit != null)
						{	model_merit = sharedmerit.get_merit();}
					}
					do
					{	run_slice();}
					while( nextHelperSlice.get() != 0 );
					if(seat == 0)
					{	
						run_mesh();
					} //single threaded
				}
				else
				{	
					seat = instance_count.decrementAndGet();
					break;
				} //user cancelled Dialog
			}
			//false decrement for very last one
			while( ( seat = instance_count.decrementAndGet() ) == 0); 
			//System.out.println(getClass().getSimpleName() + ": " + ( (seat==0)?"master":"slave" ) + 
			//					" is quitting, thread group size: " + instance_count.get() );
			clear_status();
			return;
		}
		catch( Exception e)
		{
			seat = instance_count.decrementAndGet();
			e.printStackTrace();
			System.out.println("A threading level Exception has occured, cancelling all remaining tasks, threads remaining: ~" + seat);
			nextHelperSlice.set(0);
			dialogs.clear();	
		}
		return;
	}

	void run_slice()
	{
		int my_slice = -1;
		try
		{
			int nextSlice = -1;
			do
			{
				nextSlice = nextHelperSlice.getAndSet(-1); // all others will retry every 500ms
				if(nextSlice == -1)
				{	IJ.wait(500); }
			}
			while( (nextSlice == -1) );
			if(nextSlice == 0)
			{
				nextHelperSlice.set(0);
				return;
			}
			
			
			my_slice = nextSlice; 
			boolean seat_ok = true;
			
			for(int s = 0; seat_ok && (s<taken_slices.length()); ++s)
			{ 
				//should not happen, but technically it should be ok if the same trheads restarts on the same slice
				seat_ok &= ( (s==seat) || (taken_slices.get(s) != my_slice) );
			}
			if(seat_ok)
			{	taken_slices.set(seat,my_slice); }
			//Current strategy is to wait and then switch to advance nextSlice, so dead threads would not lock the overall progress
			//However to many threads will start skipping slices
			if(do_stack)
			{
				++nextSlice;
				if(nextSlice > impDepth)
				{	nextSlice -= impDepth;}
			}
			
			if(!seat_ok) //we did not reserve a slice
			{	
				if(!waited)
				{
					System.out.println("waiting for busy slice("+my_slice+") retry in 5s");
					IJ.wait(5000); //no other thread can start a new slice as we have set it to -1
					nextHelperSlice.set(my_slice);
					waited = true;
				}
				else
				{
					System.out.println("skipping busy slice("+my_slice+") wait for 5s");
					IJ.wait(5000); //no other thread can start a new slice as we have set it to -1
					nextHelperSlice.set(nextSlice);
					waited = false;
				}
				return;
			}
						
			if(nextSlice == impSlice0)
			{ 
				completed_cycles.incrementAndGet();
			}
			//It ahould be impossible some other thread could interfer here
			//but even if he does we are most probably fine 
			
			if( (cycles > -1) && (completed_cycles.get() >= cycles) )
			{
				nextSlice = 0; //no more slices for others
			}
			waited = false;
			//System.out.println("run_slice("+my_slice+")");
			if(! very_first_slice)
			{	ensure_all_images(my_slice);}
			else
			{	very_first_slice = false;}
			
			if(reset_enhanced || (completed_cycles.get() > 0) ) 
			{	//there is no problem if someone else would set it true meanwhile
				reset_enhanced = false;				
			}
			
			if(nextHelperSlice.get() != -1) //some thread has messed with our flagged value
			{   //just cancel and quit
				taken_slices.set(seat,0);
				nextHelperSlice.set(0);
				return;
			}
			
			nextHelperSlice.set(nextSlice); //let someone else checkout the nextSlice now
			
			pull_pixels(0);
			model_merit = img_match(sim_pix,imp_pix,impWidth);
			if( (seat == 0) && (impSlice == impSlice0) )
			{	
				IJ.log("initial match of "+simulatedT+"("+ impSlice +") and "+impT+"("+ impSlice +"): " +
				model_merit + " atoms: " + numAtoms + " holes: " + numHoles);
			}
			//see if we could make use of another thread
			if( (seat+1 == instance_count.get()) && ( (seat+1) < numTasks ) &&
				(nextHelperSlice.get() != 0) && (nextHelperSlice.get() != my_slice) )
			{
				//System.out.println("thread: " + seat + " will deploy another instance");
				IJ.doCommand( getClass().getSimpleName().replaceAll("_"," ") ); 
			}	
			
			/// finally do the actual work on our slice
			IJ.resetEscape();
			if(optimize_enhanced)
			{
				show_status("thread: " + seat + " wild guesses: " + grand_runs);
				improve_model(grand_runs);
				if( optimize_Brightness )
				{	improve_Brightness(true);}	
			}
			IJ.resetEscape();
			
			
			init_mesh(mesh_from_xyz);
			if(fix_topology)
			{ fix_mesh();}
			show_mesh(mesh,enhanced);	
			
			if(optimize_fine)
			{
				show_status("thread: " + seat + " fine tuning: " + fine_runs);
				improve_positions(-1,fine_runs); //-1 .. for every Atom
				if( optimize_Brightness )
				{	improve_Brightness(true);}
				init_mesh(false);
				show_mesh(mesh,enhanced);
			}
			IJ.resetEscape();
			if(optimize_FWHM)
			{
				if( optimize_Brightness )
				{	improve_Brightness(true);}	
				improve_FWHM();
			}
			if( optimize_Brightness && !(optimize_FWHM || optimize_fine || optimize_enhanced) )
			{
				improve_Brightness(false);
			}
			
			model_merit = img_match(sim_pix,imp_pix,impWidth);
			push_pixels(0,true);
			
			
			if(seat == 0)
			{
				IJ.showProgress(0.0);
				if( (impSlice == impSlice0) )
				{	
					IJ.log("final match of "+simulatedT+"("+ impSlice +") and "+impT+"("+ impSlice +"): " +
					model_merit + " atoms: " + numAtoms + " holes: " + numHoles);
				}
				if(numTasks > 1)
				{
					IJ.log("thread slots: " + numTasks + " running: " + instance_count.get() + " busy slices: " + taken_slices);
				}
			}
			show_status("thread: " + seat + " done: " + model_merit);	
			taken_slices.set(seat,0); //Any other thread may now start on "my_slice"
			//check_floor(false);
			count_and_report_atoms(enh_pix,true);
		}
		catch(Exception e)
		{
			e.printStackTrace();
			IJ.handleException(e);
			IJ.log("malicious user interaction? cancelling all further runs");
			System.out.println("thread:" + seat + " has crashed, cancelling all remaining runs");
			nextHelperSlice.set(0);
			taken_slices.set(seat,0);
		}
			
		return; //and look for another slice
	}

	void show_status(String status_str)
	{
		if(status_str == null)	{	return;}
		Overlay po = imp.getOverlay();
		if(po != null && status_roi != null && po.contains(status_roi))
		{
			po.remove(status_roi);
		}
		status_roi = new TextRoi(10,10,status_str);
		status_roi.setPosition(impSlice);
		if(po == null)
		{
			po = new Overlay(status_roi);	
		}
		else
		{
			po.add(status_roi);
		}
		imp.setOverlay(po);
		imp.updateAndRepaintWindow();	
	}

	void clear_status()
	{
		Overlay po = imp.getOverlay();
		if(po != null && status_roi != null && po.contains(status_roi))
		{
			po.remove(status_roi);
			imp.setOverlay(po);
			imp.updateAndRepaintWindow();
		}	
	}

	void run_mesh()
	{
		//System.out.println("run_mesh()");
		master_mesh.clear();
		if(join_meshes) //only the master does the debug stuff
		{	
			//System.out.println("joining_meshes");
			for(int i = 0; i < impDepth; ++i)
			{
				final int j = ((impSlice0-1+i)%impDepth);
				mesh = meshes[j];
				mesh.incorporate(true);
				show_mesh(mesh,simulated);
			}
			ImagePlus joined = NewImage.createFloatImage(WindowManager.makeUniqueName("joined"),impWidth,impHeight,1,0);
			
			float[] pixels = (float[])joined.getStack().getPixels(1);
			Arrays.fill(pixels,Float.NaN);
			Atom[] matoms = master_mesh.fatoms;
			if(matoms != null)
			{
				for(int i = 0; i < matoms.length; ++i)
				{
					final int x = (int)(matoms[i].pos[0] + 0.5f);
					final int y = (int)(matoms[i].pos[1] + 0.5f);
					if( (x >= 0) && (x < impWidth) && (y>=0) && (y <= impHeight))
					{	pixels[x+y*impWidth] = matoms[i].observers;}
				}
			}
			
			IJ.run(joined, "Grays", "");
			show_mesh(master_mesh,joined);
			joined.show();
			IJ.log("created " + joined.getTitle());
		}
		else if(write_sig_file || write_top_file)
		{
			meshes[impSlice-1].incorporate(true);
		}
		
		if(count_atoms)
		{
			String dualname = "dual_" + enhancedT;
			if(enhancedT.equals("dual_model.tif")) //flip between the two as long as the default names are used
			{	dualname = "model.tif";}
			
			if(dual == null || !dual.getTitle().equals(dualname) )
			{	dual = WindowManager.getImage(dualname);}
			if(dual == null || !dual.isVisible() || dual.getWidth() != impWidth || dual.getHeight() != impHeight || dual.getStackSize() != impDepth)
			{
				if(dual != null && dual.isVisible())
				{	dual.close();}
				dual = NewImage.createFloatImage(dualname,impWidth,impHeight,impDepth,0);
				IJ.run(dual, "Grays", "");
			}
			
			dualSt = dual.getStack();
			dual_pix = (float[])dualSt.getPixels(impSlice);
			Arrays.fill(dual_pix,0.0f);
			if(enable_holes)
			{
				Atom[] fatoms = mesh.fatoms;
				for(int i = 0; i < fatoms.length; ++i)
				{
					float[] pos = Arrays.copyOf(fatoms[i].pos,3);
					int q = ((int)(pos[0]+0.5 + impWidth) ) % impWidth;
					int r = ((int)(pos[1]+0.5 + impHeight) ) % impHeight;
					if(hex_pixels)
					{
						int qr = imgHP.to_qr(mesh.ixyz(pos));
						q = qr%impWidth;
						r = qr/impWidth;
					}
					dual_pix[r*impWidth + q] = -0.5f * fatoms[i].intensity;
					
				}
			}
			
			Ring[] frings = mesh.frings;
			SurfaceMesh dualmesh = new SurfaceMesh(impWidth,impHeight,impSlice);
			dualmesh.hex_pixels = hex_pixels;
			dualmesh.imgHP = imgHP;
			dualmesh.periodic = periodic;
			dualmesh.inverted = !inverted; //it is supposed to be the dual
			for(int i = 0; i<frings.length; ++i)
			{
				float[] pos = Arrays.copyOf(frings[i].pos,3);
				int q = ((int)(pos[0]+0.5 + impWidth) ) % impWidth;
				int r = ((int)(pos[1]+0.5 + impHeight) ) % impHeight;
				if(hex_pixels)
				{
					int qr = imgHP.to_qr(dualmesh.ixyz(pos));
					q = qr%impWidth;
					r = qr/impWidth;
				}
				
				Atom atomi = dualmesh.put_atom(pos);
				atomi.intensity = (enable_holes?-0.5f:-1.0f) * frings[i].intensity;
				dual_pix[r*impWidth + q] = atomi.intensity;
			}
			
			dual.show();
			dualmesh.refresh();
			if(show_topology)
			{	show_mesh(dualmesh,dual);}
			dual.updateAndRepaintWindow();
			IJ.log(dualname + "(" + impSlice + ") with " + frings.length + " atoms" + (enable_holes?(" and " + mesh.fatoms.length + " holes"):""));
			
		}
		
		
		if(write_sig_file)
		{	write_sig(master_mesh);}
		
		if(write_top_file)
		{	write_top();}
		
	}



	void improve_model(int full_trials)
	{
		boolean keep_going = false;
		boolean single_atom = false;
		if(full_trials < 0)
		{	
			full_trials = 5000;
			if(seat == 0)
			{	IJ.log("looking for one more atom " + (enable_holes?"/hole":""));}	
			keep_going = true;
			single_atom = true;
		}
		
		
		final boolean[] visited = new boolean[impArea];
		if(hex_pixels)
		{
			for(int i = 0; i < impArea; ++i)
			{	visited[i] = !imgHP.is_inside(i);}
		}
		else
		{
			//for(int i = 0; i < impArea; ++i)
			//{	visited[i] = Float.isNaN(imp_pix[i]);}
		}
		int visited_pts = 0;
		int visited_resets = 0;
		double tentative_merit = 0;
		do
		{
			if(seat == 0)
			{	IJ.showStatus("commencing " + full_trials + " grand trials, hit Esc to skip");}
			
			int accepted = 0;
			int rejected = 0;
			int stable = 0;
			int actual_trials = 0;
			
			for(int ft = 0; ft < full_trials; ++ft)
			{
				int trial = 0;
				++actual_trials;
				
				int nextpos = -1;
				double nextweight = Double.NaN;
					
				do
				{
					if(visited_pts >= (hex_pixels ? imgHP.hpMax : validArea))
					{
						if(hex_pixels)
						{
							//Arrays.fill(visited,true);
							for(int i = 0; i < impArea; ++i)
							{	visited[i] = !imgHP.is_inside(i);}
						}
						else
						{
							Arrays.fill(visited,false);
						}
						visited_pts = 0;
						++visited_resets;		
					}
					
					//check out one random point before picking the max once
					//if we are looking for one more atom, always go with the max
					if(single_atom || trial==1)
					{
						nextpos = get_max_pos(dev_pix, visited);
						if(nextpos < 0)
						{ 
							///This can very rarely happen on noisy images
							//System.out.println("Could not find any unvisted maxima > 0");
							++trial;
							continue;
						}
						if(hex_pixels && !imgHP.is_inside(nextpos))
						{	System.out.println("nextpos from get_max was not inside hex domain");}	
					}
					else
					{
						nextpos = -1;
						int nextpick = random.nextInt((hex_pixels ? imgHP.hpMax : validArea)-visited_pts);
						final int nextpick0 = nextpick;
						do
						{
							if(!Float.isNaN(imp_pix[++nextpos]))
							{
								nextpick -= ( visited[nextpos]?0:1 );
							}
						}
						while(nextpick >= 0);
						if(nextpos < 0)
						{
							System.out.println(	"Random pick failed, nextpos: " + nextpos +
												" nextpick: " + nextpick + " nextpick0: " + nextpick0 + " visited_pts: " + visited_pts);			
							++trial;
							continue;
						}
						/*
						if(hex_pixels && !imgHP.is_inside(nextpos))
						{	System.out.println("nextpos from random pick was not inside hex domain");}
						if(!hex_pixels && Float.isNaN(imp_pix[nextpos]))
						{	System.out.println("nextpos from random pick was on NaN area");}
						*/ 
					}
					
					
					nextweight = estimate_weight(enh_pix,nextpos);
					//mark the entire solid crosssection of the tentative atom as visited
					//if(nextweight != -1.0)
					{
						final int px = nextpos % impWidth;
						final int py = nextpos / impWidth;
						
						final int rmax = (int)Math.ceil(0.5*atomSolidD);
						final double wsqr = 0.25*atomSolidD*atomSolidD;
						for(int y = -rmax; y<=rmax; ++y)
						{
							for(int x = -rmax; x<=rmax; ++x)
							{
								
								if(!hex_pixels)
								{
									if(x*x+y*y > wsqr)
									{	continue;}
								} 
								else //hex_pixels is true
								{
									if(imgHP.peucd2qr(px,py,px+x,py+y)>wsqr)
									{continue;}
								}
								
								int ax = px+x;
								int ay = py+y;
								if(!hex_pixels) 
								{
									if( periodic)
									{
										ax = (ax+impWidth)%impWidth;
										ay = (ay+impHeight)%impHeight;
									}
									else if( (ax < 0) || (ay < 0) ||
											 (ax >= impWidth) || (ay >= impHeight) )
									{	continue;}
								}
								else //hex_pixels = true
								{
									final int pqr = imgHP.periodic_qr_safe(ax,ay);
									ax = pqr%impWidth;
									ay = pqr/impWidth;
								}
								
								
								final int pos = ax+ay*impWidth;
								if(!visited[pos])
								{
									visited[pos] = true; //dont get stuck on hot pixels
									++visited_pts;
								}
								
							}
						}
					}
					/*else
					{
						visited[nextpos] = true; //dont get stuck on hot pixels
						++visited_pts;
					}*/
					
				}
				while( (++trial < 20) && Double.isNaN(nextweight) );
				
				if( Double.isNaN(nextweight) || (nextpos == -1) ) 
				{
					++stable;
					continue;
				}
				
				if( (!enable_holes) && (nextweight < 0.0))
				{
					throw new RuntimeException("estimated weight: " + nextweight + " although enable_holes: " + enable_holes);
				}
				
				
				if(!optimize_weights)
				{	
					if(nextweight * avg_peak_weight > 0.0) 
					{	nextweight = avg_peak_weight;}
					else if (nextweight * avg_hole_weight > 0.0)
					{	nextweight = avg_hole_weight;}	
				}
				
				
				double repelled = put_down_atom(enh_pix,impWidth,nextpos,nextweight);
				
				if(nextweight != 0.0) //we did actually place an atom and let improve_positions take care of everything
				{
					int newpos = nextpos;
					int imppos = nextpos;
					int moves = 0; //we might be jumping forth and back between 2 equal spots
					int onA = numAtoms;
					do
					{
						newpos = imppos;
						//impos may become -1 if the atom has crashed
						/*
						if(enh_pix[newpos] == 0.0f )
						{
							throw new RuntimeException("improve_positions has previously fucked up an atom in move: " + moves);
						}
						*/ 
						imppos = improve_positions(newpos,2); //the second pass will compensate for kicked out atoms
						/*
						if(enh_pix[imppos] == 0.0f)
						{
							throw new RuntimeException("improve_positions has just fucked up an atom in move: " + moves);
						}
						if(imppos != newpos && (enh_pix[newpos] != 0.0f) )
						{
							throw new RuntimeException("improve_positions has just duplicated an atom in move: " + moves);
						}
						*/ 
					}
					//limit steps as well as number of kicked out atoms
					while(imppos != newpos && (++moves <= 3) && (imppos >= 0) && (numAtoms == onA) ); 
				}
				else //we have to do the simulation and calculate the merit
				{
					sim_pix = sparse_floor_fold(enh_pix,impWidth, modelFloor, a_pix,atomWidth);
					model_merit = img_match(sim_pix,imp_pix,impWidth);
				}
				
				push_pixels(0,true);
				++accepted;
				
				if(seat == 0)
				{
					IJ.showProgress(ft,full_trials);
				}
				if(IJ.escapePressed() || single_atom)
				{ 
					keep_going = false;
					break;
				}
				/*
				if(repelled != 0.0 && nextweight != 0.0)
				{
					System.out.println("just moved an atom to x,y: " + (nextpos%impWidth) + "," + (nextpos/impWidth));
					System.out.println("new weight: " + nextweight);
					System.out.println("last weight: " + repelled);
					keep_going = false;
					break;
				}
				*/
			
			}
			keep_going &= (accepted > 0);
			final int stableratio = (actual_trials > 0) ? 100*stable/actual_trials : 0;
			final int visitedratio = 100*visited_resets + ((visited_pts > 0) ? 100*visited_pts/(hex_pixels ? imgHP.hpMax : validArea) : 0);
			if( (seat == 0) && (!do_stack) )
			{	
				IJ.log("thread: " + seat + " finished " + actual_trials + " grand trials stable: " + stableratio +
				 "% visited: " + visitedratio + "% slice("+ impSlice +") merit: " + model_merit + " atoms: " + numAtoms + " holes: " + numHoles);
			}
			
		
		}
		while(keep_going);
	}
	
	void init_mesh(boolean from_xyz)
	{
		mesh = meshes[impSlice-1];
		
		int count = 0;
		if(from_xyz)
		{
			mesh.clear();
			count = load_xyz();
			mesh.refresh();
			if( count == 0 )
			{
				count = 0;
				for(int i = 0; i < impArea; ++i)
				{
					if(enh_pix[i] > 0)
					{	++count;}
				
				}
			}
		
		}
		
		if(count == 0)
		{
			mesh.clear();
			for(int i = 0; i < impArea; ++i)
			{
				final float eval = enh_pix[i];
				if(eval>0.0f)
				{
					++count;
					float[] pos = null;
					if(!hex_pixels)
					{
						final float x = i%impWidth;
						final float y = i/impWidth;
						pos = new float[]{x,y,0.0f};  
					}
					else
					{
						final int[] xyz = imgHP.to_xyz(i);
						pos = new float[]{xyz[0],xyz[1],xyz[2]};
					}
					Atom ma = mesh.put_atom( pos );
					ma.intensity = eval;
				}
			}
		}
		numAtoms = count;
		mesh.refresh();
	}
	
	int load_xyz()
	{
		int count = 0;
		FileInfo fi = imp.getOriginalFileInfo();
		String path = (fi!=null)?fi.directory:IJ.getDirectory("current");
		if(path == null) path = IJ.getDirectory("current");//just in case imp did not come from disk
		
		String fname = imp.getShortTitle() + "#" + impSlice  + ".top";
		String xyzT = null;
		File xyzF = new File(path+fname);
		if(! (xyzF.isFile() && xyzF.canRead()))
		{
			OpenDialog od = new OpenDialog("Open xyz file for initial mesh in (" + impSlice + ")",
			path,fname);
			xyzT = od.getPath();
		}
		if(xyzT == null)
		{	return 0;} //just go ahead with scanning the model image
		try
		{
			
			BufferedReader reader = new BufferedReader(new FileReader(xyzT));
			String line;
			reader.readLine(); //number of atoms
			reader.readLine(); //mandatory comment line
			while (reader.ready())
			{
				line = reader.readLine();
				if(line == null)
				{	break;}
				String[] words = line.split("\\s+");
				if(words.length >= 4)
				{
					String tag = words[0]; //Chemical Element or generic "ATOM"
					if(tag.equals("ATOM") || tag.equals("C"))
					{
						++count;
						float x = Float.valueOf(words[1]);
						float y = Float.valueOf(words[2]);
						float z = Float.valueOf(words[3]);
						float[] pos = new float[]{x,y,z};
						Atom ma = mesh.put_atom( pos );
						ma.intensity = 1.0f;
						enh_pix[(int)(y+0.5)*impWidth+(int)(x+0.5)]=1.0f;  
					}
					
				}
			  
			  
			}
			reader.close();
			IJ.log("imported " + count + " atoms from " + xyzT);
			return count;
		}
		catch (Exception e)
		{
			System.out.println("failed to read: " + xyzT); 
			e.printStackTrace();
			return 0;
		}
	}
	
	
	void show_mesh(SurfaceMesh m, ImagePlus b)
	{
		m.drawMarks(b, show_topology, show_topology, show_topology, show_topology);
		if(b == enhanced)
		{
			RoiManager roim = RoiManager.getRawInstance();
			//spare the effort if the user is not using the RoiManager anyways
			//Roi[] rois = b.getOverlay().toArray();
			//if( (rois != null) && (rois.length>0) )
			//only push actually drawn rois to RoiManager
			if((roim != null) &&  (roim.isShowing()) && (b.getOverlay().size() > 0) )
			{
				
				IJ.run(b,"To ROI Manager", "");
				roim.runCommand(b,"show all without labels");
				
				/* //This makes RoiManger Window Flicker
				roim.reset();
				for(int i = 0; i < rois.length; ++i)
				{	roim.addRoi(rois[i]);}
				*/ 
			}	
		}	
	}
	
	void fix_mesh()
	{
		//do only one action at a time, the user has to decide on multiple runs
		Atom[] fatoms = mesh.fatoms;
		Bond[] fbonds = mesh.fbonds;
		Ring[] frings = mesh.frings;
		boolean action_taken = false;
		
		if(!inverted) // the typical complicated use case
		{
			//System.out.println("scanning bonds: " + fbonds.length);
			//1 scan for suspicious bonds
			for(int i = 0; i < fbonds.length; ++i)
			{
				Bond bondi = fbonds[i];
				if(bondi.detached)
				{	continue;}
				if(!bondi.is_interior)
				{	continue;}
				if( (bondi.left_ring == null) || (bondi.right_ring == null) )
				{	continue;}
				if((bondi.left_ring.edges.size() != 5) || (bondi.right_ring.edges.size() != 5))
				{	continue;}
				//System.out.print(".");
				Atom a1 = bondi.a1;
				Atom a2 = bondi.a2;
				//check if a bond has cut a hexagon via two interstitial atoms
				if((a1.rings.size() == 3) && (a2.rings.size() == 3))
				{	
					int a1r7 = 0;
					int a2r7 = 0;
					for(int j = 0; j < 3; ++j)
					{
						if(a1.rings.get(j).edges.size()>=7)
						{	++a1r7;}
						if(a2.rings.get(j).edges.size()>=7)
						{	++a2r7;}
					}
					if( (a1r7!=1) || (a2r7!=1) )
					{	continue;}
					enh_pix[a1.ind] = 0.0f;
					enh_pix[a2.ind] = 0.0f;
					mesh.remove_atom(a1);
					mesh.remove_atom(a2);
					action_taken = true;
					System.out.println("removed 7557 bond near x,y : " + (int)bondi.pos[0] + "," + (int)bondi.pos[1]);
					break;
				}
				//scan for a Z of 5 rings
				else if((a1.rings.size() == 4) && (a2.rings.size() == 4) )
				{
					Ring r1 = null;
					Ring r2 = null;
					for(int k = 0; (k < 4) && ( (r1==null) || (r2==null) ); ++k)
					{
						if( (r1==null) && (!a1.rings.get(k).has_edge(bondi)) )
						{	r1 = a1.rings.get(k);}
						if( (r2==null) && (!a2.rings.get(k).has_edge(bondi)) )
						{	r2 = a2.rings.get(k);}
					}
					//if any of them stays null the topology is no longer 2D! 
					if( (r1.edges.size()<=5) && (r2.edges.size()<=5) )
					{
						// We have Z shape of 4 5-rings, there are two atoms missing
						// replace a1 and a2 by the pair of midpoint positions along
						// the edges of 5-rings at the ends
						System.out.println("fixing (4)5-ring Z-shape around x,y: " + (int)bondi.pos[0] + "," + (int)bondi.pos[1]);
						for(int k=0; k < 4; ++k)
						{
							if(r1.has_edge(a1.edges.get(k)))
							{
								//Ring ro = a1.edges.get(k).get_other_ring(r1);
								Bond bk = a1.edges.get(k);
								//if(ro != null)
								{
									float mx = bk.pos[0];//0.5f*(r1.pos[0] + a1.pos[0]);
									float my = bk.pos[1];//0.5f*(r1.pos[1] + a1.pos[1]);
									float mz = bk.pos[2];//0.5f*(r1.pos[2] + a1.pos[2]);
									float[] mp = {mx,my,mz};
									Atom ma = mesh.put_atom(mp);
									ma.intensity = a1.intensity;
									enh_pix[ma.ind] = ma.intensity;
								}
							}
							if(r2.has_edge(a2.edges.get(k)))
							{
								//Ring ro = a2.edges.get(k).get_other_ring(r2);
								Bond bk = a1.edges.get(k);
								//if(ro != null)
								{
									float mx = bk.pos[0];//0.5f*(r2.pos[0] + a2.pos[0]);
									float my = bk.pos[1];//0.5f*(r2.pos[1] + a2.pos[1]);
									float mz = bk.pos[2];//0.5f*(r2.pos[2] + a2.pos[2]);
									float[] mp = {mx,my,mz};
									Atom ma = mesh.put_atom(mp);
									ma.intensity = a2.intensity;
									enh_pix[ma.ind] = ma.intensity;
								}
							}
						}
						enh_pix[a1.ind] = 0.0f;
						mesh.remove_atom(a1);
						enh_pix[a2.ind] = 0.0f;
						mesh.remove_atom(a2);
						action_taken = true;
						break;
					}
					
				} 
			}
			mesh.refresh();
			
			frings = mesh.frings;
			fbonds = mesh.fbonds;
			fatoms = mesh.fatoms;
			//only proceed with atoms if the bonds checked out
			if(!action_taken)
			{
				//System.out.println("scanning atoms: " +fatoms.length);
				for(int i = 0; i < fatoms.length; ++i)
				{
					Atom atomi = fatoms[i];
					if(atomi.detached)
					{ continue;}
					//System.out.print(".");
					if(atomi.rings.size() == 4)
					{
						int r4 = 0;
						int r5 = 0;
						int r6 = 0;
						int r7 = 0;
						int r8 = 0;
						int rx = 0;
						for(int j = 0; j < atomi.rings.size(); ++j)
						{
							Ring ringj = atomi.rings.get(j);
							if(ringj.edges.size() == 4)
							{	++r4;}
							else if(ringj.edges.size() == 5)
							{	++r5;}
							else if(ringj.edges.size() == 6)
							{	++r6;}
							else if(ringj.edges.size() == 7)
							{	++r7;}
							else if(ringj.edges.size() == 8)
							{	++r8;}
						}
						rx = atomi.rings.size()-(r4+r5+r6+r7+r8);
						//split the atom, and push fragments into non 5 rings
						//except into 6 rings if they are the only local 6 ring
						if( (r5+r4>=2) && (r7+r8 != 2) && (r4 != 2) )
						{
							final float f = (6-(r5+r4+((r6==1)?1:0)))/6.0f;
							System.out.println("seeding atoms around 4-fold atom at x,y: " + (int)atomi.pos[0] + "," + (int)atomi.pos[1]);
							int created = 0;
							for(int j = 0; j < atomi.rings.size(); ++j)
							{
								final Ring ringj = atomi.rings.get(j);
								final int rjes = ringj.edges.size();
								if(ringj.edges.size() != 5 //&& (! ( (rjes==6) && (r6 == 1) ) )
								  )
								{
									float mx = (f*atomi.pos[0] + (1.0f-f)*ringj.pos[0]);
									float my = (f*atomi.pos[1] + (1.0f-f)*ringj.pos[1]);
									float mz = (f*atomi.pos[2] + (1.0f-f)*ringj.pos[2]);
									float[] mp = {mx,my,mz};
									Atom ma = mesh.put_atom(mp);
									ma.intensity = atomi.intensity;
									enh_pix[ma.ind] = ma.intensity;
									++created;
									System.out.println("created atom at x,y: " + (int)ma.pos[0] + "," + (int)ma.pos[1]);
								}
							}
							if(created > 1)
							{
								enh_pix[atomi.ind] = 0.0F;
								mesh.remove_atom(atomi);
								System.out.println("remove central atom at x,y: " + (int)atomi.pos[0] + "," + (int)atomi.pos[1]);
							}
							if(created > 0 )
							{	break;}
						}		
					}
					else if(atomi.edges.size()==1 && atomi.rings.isEmpty())
					{
						float[] ip = atomi.pos;
						Atom other = atomi.neighbors.get(0);
						float[] op = other.pos;
						if(!other.rings.isEmpty())
						{
							for(int r=0; r < other.rings.size(); ++r)
							{
								float[] rp = other.rings.get(r).pos;
								if(mesh.distsqr(rp,op) > mesh.distsqr(rp,ip))
								{
									enh_pix[atomi.ind] = 0.0f;
									mesh.remove_atom(atomi);
									System.out.println("removed enclosed dangling atom at x,y: " + (int)atomi.pos[0] + "," + (int)atomi.pos[1]);
									r=other.rings.size();
									break;
								}
							}
						}
						
					}
					/*
					else if(atomi.rings.size() == 3)
					{
						int r5 = 0;
						int r7 = 0;
						for(int j = 0; j < atomi.rings.size(); ++j)
						{
							Ring ringj = atomi.rings.get(j);
							if(ringj.edges.size() == 5)
							{	++r5;}
							if(ringj.edges.size() == 7)
							{	++r7;}
						}
						if( (r5 == 2) && (r7 == 1) )
						{
							enh_pix[atomi.ind] = 0.0F;
							mesh.remove_atom(atomi);
							IJ.log("removed 557 atom at x,y : " + (int)atomi.pos[0] + "," + (int)atomi.pos[1]);
						}	
					}
					else if( (atomi.edges.size() == 2) && (atomi.rings.size()==0) )
					{
						if( (!atomi.neighbors.get(0).detached) && (!atomi.neighbors.get(1).detached) )
						{ //kill every second atom in freestanding line
							enh_pix[atomi.ind] = 0.0F;
							mesh.remove_atom(atomi);
							IJ.log("removed free linear atom at x,y : " + (int)atomi.pos[0] + "," + (int)atomi.pos[1]);
							continue;
						}
					}
					*/ 
					else if( (atomi.edges.size() == 2) && (atomi.rings.size() == 2) )
					{
						
						
						
						
						final float[] v0 = atomi.vectorTo(atomi.neighbors.get(0));
						final float[] v1 = atomi.vectorTo(atomi.neighbors.get(1));
						final float sk = v0[0]*v1[0]+v0[1]*v1[1]+v0[2]*v1[2];
						final float sk0 =  v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2];
						final float sk1 =  v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
						final double c = (sk/Math.sqrt(sk0*sk1));
						if(c < -0.85)
						{
							enh_pix[atomi.ind] = 0.0f;
							mesh.remove_atom(atomi);
							System.out.println("removed nearly straight linear atom at x,y : " + (int)atomi.pos[0] + "," + (int)atomi.pos[1]);
							break;
						}	
						else //if( atomi.rings.size() == 2)
						{
							
							//System.out.println("keeping kinked linear atom at x,y : " + (int)atomi.pos[0] + "," + (int)atomi.pos[1]);
							Ring r0 = atomi.rings.get(0);
							Ring r1 = atomi.rings.get(1);
							if(r0.edges.size() == r1.edges.size())
							{ //most likely an extraneous atom 
								enh_pix[atomi.ind] = 0.0f;
								mesh.remove_atom(atomi);
								System.out.println("removed extraneous atom at x,y : " + (int)atomi.pos[0] + "," + (int)atomi.pos[1]);
								break;
							}
							
							
							Ring rs = r0;
							Ring rb = r1;
							if(r0.edges.size()>r1.edges.size())
							{
								rb = r0; //bigger ring
								rs = r1; //smalller ring
							}
							
							//place an atom inside the bigger ring
							///This needs to become smarter before it does more good than bad
							/*
							if(rb.edges.size() > rs.edges.size()+1)
							{
								float mx = 2*rb.pos[0] - atomi.pos[0];
								float my = 2*rb.pos[1] - atomi.pos[1];
								float mz = 2*rb.pos[2] - atomi.pos[2];
								float[] mp = {mx,my,mz};
								Atom ma = mesh.put_atom(mp);
								
								ma.intensity = atomi.intensity;
								if( (ma.ind >= 0) && (ma.ind < enh_pix.length) )
								{	
									enh_pix[ma.ind] = ma.intensity;
									System.out.println("seeded atom at x,y : " + (int)ma.pos[0] + "," + (int)ma.pos[1]);
								}
								else
								{	
									System.out.println("Failed to seeded atom at x,y : " + (int)ma.pos[0] + "," + (int)ma.pos[1]);
									mesh.remove_atom(ma);
								}
								break;
							}
							*/ 
							
						
						}
					
					}				
				}
				
				mesh.refresh();
				frings = mesh.frings;
				fbonds = mesh.fbonds;
				fatoms = mesh.fatoms;
			}
			/* hmm maybe later we do 
			boolean any_fix = false;
			do
			{		
				any_fix = false;
				mesh.refresh();
				frings = mesh.frings;
				fbonds = mesh.fbonds;
				fatoms = mesh.fatoms;
				for(int i=0; i < frings.size; ++i)
				{
					Ring ri = frings[i];
					if(ri.edges.size()!=6)
					{
						
						
					}	
				
				
				
				}
			
			}
			while(any_fix);
			*/
		}
		else //inverted == true
		{
			for(int i=0; i < frings.length; ++i)
			{
				Ring ri = frings[i];
				if(ri.edges.size()>4) //no "actual" atom could have more than 4 bonds
				{
					float[] mp = Arrays.copyOf(ri.pos,3);
					Atom ma = mesh.put_atom(mp);
					ma.intensity = ri.vertices.get(0).intensity;
					enh_pix[ma.ind] = ma.intensity;	
				}	
			}
			
		}
	
		mesh.refresh();
		fatoms = mesh.fatoms;
		numAtoms = fatoms.length;
		sim_pix = sparse_floor_fold(enh_pix,impWidth, modelFloor, a_pix,atomWidth);
		push_pixels(0,true);
		
		
	}
	
	
	
	
	
	
	void improve_FWHM()
	{
		if(numAtoms == 0)
		{	return;}
		
		final int area = atomWidth*atomHeight;
		
		double best_atomFWHM = atomFWHM;
		double best_merit = model_merit;
		double f = 0.9;
		double fmax = 0.999;
		boolean shrinking = random.nextBoolean();
		if(seat == 0)
		{
			IJ.showStatus("FWHM: " + atomFWHM);
			IJ.showProgress(model_merit);
			IJ.log("optimizing FWHM: " + atomFWHM + " match: " + model_merit);
		}
		while(f < fmax)
		{
			atomFWHM = best_atomFWHM*(shrinking?f:1/f);	
			init_atom(a_pix);
			sim_pix = sparse_floor_fold(enh_pix,impWidth, modelFloor,a_pix,atomWidth);
			double tentative_merit = img_match(sim_pix,imp_pix,impWidth);
			
			if(tentative_merit > best_merit)
			{
				best_merit = tentative_merit;
				best_atomFWHM = atomFWHM;
				if(seat == 0)
				{	
					IJ.showStatus("FWHM: " + atomFWHM);
					IJ.showProgress(best_merit);
				}
				push_pixels(0,true);
				//never accelarate as we expect to be relatively close anyways	
			}
			else
			{
				atomFWHM = best_atomFWHM;
				shrinking = !shrinking;
				f = Math.sqrt(f);//carefull slow down
			}
			if(IJ.escapePressed())
			{ 
				break;
			}
		}
		atomFWHM = best_atomFWHM;
		init_atom(a_pix);
		
		sim_pix = sparse_floor_fold(enh_pix,impWidth, modelFloor,a_pix,atomWidth);
		model_merit = img_match(sim_pix,imp_pix,impWidth);
		push_pixels(0,true);
		atom.updateAndRepaintWindow();
		if(seat == 0)
		{	IJ.log("new FWHM: " + atomFWHM + " new match: " + model_merit);}	
	}
	
	void measure_FWHM(float[] atom)
	{
		final int area = atomWidth*atomHeight;
		float[] tmp_atom = new float[area];
		double best_atomFWHM = atomFWHM;
		init_atom(tmp_atom);
		double best_match = img_match(tmp_atom,atom,atomWidth);;
		double f = Math.sqrt((1 + Math.abs(best_match) )/2);
		double fmax = 0.999;
		boolean shrinking = random.nextBoolean();
		//if(seat == 0)
		//{	IJ.log("measuring FWHM in " + atomT + " starting value: " + atomFWHM + " match: " + best_match);}
		boolean any_change = false;
		while(f < fmax)
		{
			atomFWHM = best_atomFWHM*(shrinking?f:1/f);	
			init_atom(tmp_atom);
			double tentative_match = img_match(tmp_atom,atom,atomWidth);
			
			if(tentative_match > best_match)
			{
				any_change = true;
				best_match = tentative_match;
				best_atomFWHM = atomFWHM;
				if(seat == 0)
				{
					IJ.showStatus("FWHM: " + atomFWHM);
					IJ.showProgress(best_match);
				}
			}
			else
			{
				atomFWHM = best_atomFWHM;
				shrinking = !shrinking;
				f = Math.sqrt(f);//carefull slow down
			}
			if(IJ.escapePressed())
			{ 
				break;
			}
		}
		if(any_change)
		{
			atomFWHM = best_atomFWHM;
			//if(seat == 0)
			//{	IJ.log("new FWHM: " + best_atomFWHM + " new match: " + best_match);}
		}	
	}
	
	void improve_Brightness(boolean silent)
	{
		if(numAtoms + numHoles < 1)
		{	return;}
		final int area = hex_pixels ? imgHP.hpMax : impArea;
		double f = 0.9;
		double fmax = 0.99;
		boolean fok = (f < fmax);
		double model_var = get_var(dev_pix,0.0f,imgHP);
		if(model_var == 0.0)
		{	return;}
		double model_std = Math.sqrt(model_var);
		if(  (!silent) && (seat == 0) ) //
 		{	IJ.log("Adjusting global atom brightness, current std in " + deviationsT + "(" + impSlice + "): " + model_std);}
		boolean shrinking = random.nextBoolean();
		double scale = shrinking?f:1.0/f;
		
		apply_scale((float)scale);
		while(fok)
		{
			sim_pix = sparse_floor_fold(enh_pix,impWidth, modelFloor,a_pix,atomWidth);
			for(int ind = 0; ind < area; ++ind)
			{
				final int i = hex_pixels?imgHP.hp[ind]:ind;
				final float ival = imp_pix[i];
				//dev_pix[i] = Float.isNaN(ival)?0.0:ival-sim_pix[i];
				dev_pix[i] = ival-sim_pix[i];
			}
			final double tentative_var = get_var(dev_pix,0.0f,imgHP);		
			if(tentative_var >= model_var)
			{
				shrinking = !shrinking;
				double correct_last = 1.0/scale;
				f=Math.sqrt(f);
				scale = shrinking?f:1.0/f;
				if( (fok = f < fmax) ) 
				{
					correct_last *= scale;
				}
				apply_scale((float)correct_last); //and also apply the new scale if there will be another iteration
			}
			else //just keep going
			{
				model_var = tentative_var;
				apply_scale((float)scale);
			} 		
		}
		push_pixels(0,true);
		model_std = Math.sqrt(model_var);
		
		
		if(  (!silent) && (seat == 0) ) //
		{
			IJ.log("new std in " + deviationsT + "(" + impSlice + "): " + model_std);
		}
	}
	
	void apply_scale(final float scale)
	{
		double enhnorm = 0.0;
		final int area = hex_pixels ? imgHP.hpMax : impArea;
		for(int ind=0; ind < area; ++ind)
		{
			final int i = hex_pixels?imgHP.hp[ind]:ind;
			final float enval = enh_pix[i];
			if(enval>0.0f || enval < 0.0f)
			{
				final float nenval = enval*scale;
				enh_pix[i] = nenval;
				enhnorm += nenval;
			} 
		}
		avg_peak_weight *= scale;
		avg_hole_weight *= scale;
		
		modelFloor = ((impTotal - enhnorm)/(hex_pixels ? imgHP.hpMax : validArea));
		//check_floor(false);
	}
	
	double get_var(float[] any_pix, float avg, HexPixels anyHP)
	{
		if(Float.isNaN(avg)) return Double.NaN;
		double sum2 = 0.0;
		final int area = any_pix.length;
		int varea = 0;
		if(anyHP == null)
		{
			for(int ind=0; ind < area; ++ind)
			{
				final float v = any_pix[ind]-avg; 
				if(!Float.isNaN(v))
				{
					sum2 += v*v;
					++varea;
				}	
			}
			if(varea == 0) return Double.NaN;	
		}
		else
		{
			for(int h=0; h < anyHP.hpMax; ++h)
			{
				final float v = any_pix[anyHP.hp[h]]-avg; 
				sum2 += v*v;
				
			}
			
		}
		
		return (sum2/(anyHP==null?varea:anyHP.hpMax));
	}
	
	double get_avg(float[] any_pix, HexPixels anyHP)
	{
		double sum = 0.0;
		final int area = any_pix.length;
		int varea = 0;
		if(anyHP == null)
		{
			if(hex_pixels)
			{	System.out.println("missing hexpixels for calculating average in hex_mode");}
			for(int ind=0; ind < area; ++ind)
			{
				final float v = any_pix[ind]; 
				if(!Float.isNaN(v))
				{
					sum += v;
					++varea;
				}
			}
			validArea = varea;
			if(varea == 0) return Double.NaN;	
		}
		else
		{
			
			for(int h=0; h < anyHP.hpMax; ++h)
			{
				final float v = any_pix[anyHP.hp[h]]; 
				sum += v;	
			}
			
		}
		return (sum/(anyHP==null?varea:anyHP.hpMax));
	}
	
	int get_max_pos(float[] any_pix, boolean[] visited)
	{
		int max_pos = -1;
		float max = 0.0f;
		final int area = any_pix.length;
		for(int ind = 0; ind < area; ++ind)
		{
			if(!visited[ind])
			{
				final float v = (enable_holes?Math.abs(any_pix[ind]):any_pix[ind]); 
				if( v > max )
				{
					max = v;
					max_pos = ind;
				}
			}
		}
		return (max_pos);
	}
	
	
	
	int improve_positions(final int ind0, int full_runs)
	{
		final boolean single_atom = (ind0 > -1);
		final int oldnumAtoms = numAtoms;
		if(!single_atom) 
		{
			/* //ok simply treast all atoms carefully if D < 4
			if(atomSolidD < 4)
			{
				IJ.log("thread: "+seat+" Warning, criticaly small solid diameter for nudging : " + atomSolidD);
			}
			*/ 
			if (numAtoms + numHoles < 1 )
			{
				IJ.log("thread: "+seat+" Sorry, no atoms no fine tuning");
				return -1; //That will crash improve_model
			}
		}
		
		final double atomSolidD_backup = atomSolidD;
		final double atomSolidD_tiny = 1;//(1 + atomSolidD)/2; //make colissions unlikely//atoms can wander inside the death zones
		boolean cancelled = false;
		boolean endless = false;
		if(full_runs < 0)
		{
			full_runs = 10;
			endless = true;
		}
		
		int nA = 0;
		int[] positions = null;
		float[] weights = null;
		
		do
		{
			if(!single_atom)
			{
				if(nA != numAtoms+numHoles)
				{
					nA = numAtoms+numHoles;
					positions = new int[nA];
					weights = new float[nA];
					int p = 0;
					for(int ind = 0; ind < impArea; ++ind)
					{
						float eval = enh_pix[ind];
						if( (eval > 0.0f) || (eval < 0.0f) )
						{
							positions[p]=ind;
							weights[p]=eval;
							++p;
						}
					}
				}
			}
			else if(nA==0)//single_atom == true
			{
				nA = 1;
				positions = new int[nA];
				weights = new float[nA];
				positions[0] = ind0;
				weights[0] = enh_pix[ind0];
			}
			
			
			if( (!single_atom) && (seat == 0) )
			{
				IJ.showStatus("commencing " + full_runs + " passes of fine tuning, hit Esc to skip");
			}
			
			for(int fr = 0; (fr < full_runs) && (!cancelled) ; ++fr)
			{
				atomSolidD = atomSolidD_tiny; //shrink atoms and "set up traps" 
				if(!single_atom)
				{
					for(int i = 1; i < nA; ++i)
					{
						final int j = random.nextInt(i+1);
						if(i!=j)
						{
							final int p = positions[j];
							positions[j] = positions[i];
							positions[i] = p;
							final float w = weights[j];
							weights[j] = weights[i];
							weights[i] = w;	
						}
					}
				}
				
				
				//step1 nudge positions within wobble radius
				
				for(int a = 0; (a < nA) && (!cancelled) ; ++a)
				{	
					int posa = positions[a];
					float weighta = weights[a];
					//some other locally optimal atom has crashed into this one
					if( (posa == -1) || (enh_pix[posa] == 0.0f) ) 
					{
						positions[a] = -1;
						continue;
					}
					
					
					final int px = posa%impWidth;
					final int py = posa/impWidth;
					
					enforce_atom(enh_pix,impWidth,posa,0.0); 
					double weight0 = single_atom?weighta:estimate_weight(enh_pix, posa);
					final int rmax = (int)Math.ceil(wobble);
					final double wsqr = wobble*wobble;
					
					
					for(int ry = 0; ry<=2*rmax; ++ry)
					{
						final int y = ((ry+1)/2)*((ry%2)==0?1:-1);
						for(int rx = 0; rx<=2*rmax; ++rx)
						{
							final int x = ((rx+1)/2)*((rx%2)==0?1:-1);
							if(!hex_pixels)
							{
								if(x*x+y*y > wsqr)
								{	continue;}
							} 
							else //hex_pixels is true
							{
								if(imgHP.peucd2qr(px,py,px+x,py+y)>wsqr)
								{continue;}
							}
							
							int ax = px+x;
							int ay = py+y;
							if(!hex_pixels) 
							{
								if( periodic)
								{
									ax = (ax+impWidth)%impWidth;
									ay = (ay+impHeight)%impHeight;
								}
								else if( (ax < 0) || (ay < 0) ||
										 (ax >= impWidth) || (ay >= impHeight) )
								{	continue;}
							}
							else //hex_pixels = true
							{
								final int pqr = imgHP.periodic_qr_safe(ax,ay);
								ax = pqr%impWidth;
								ay = pqr/impWidth;
							}
							
							
							final int pos = ax+ay*impWidth;
							
							//dont kick out another atom just for checking out
							//dont even try to leave the actual data
							if( (enh_pix[pos] != 0.0f) || (Float.isNaN(imp_pix[pos])))
							{	continue;} 
							
							final double nweight = estimate_weight(enh_pix, pos);
							if( (nweight*weight0>=0.0) && (Math.abs(nweight) > Math.abs(weight0)) )
							{
								weight0 = nweight;
								positions[a] = pos;
								//set the estimated weights only for fresh atoms
								//dont interfer but only preced the global weight optimization 
								if(single_atom && (!optimize_weights) )
								{	weights[a] = (float)nweight;}
							}
							
							
							
						}
					}
					/*
					if( posa != positions[a] && (enh_pix[posa] != 0.0f) )
					{
						enforce_atom(enh_pix,impWidth,posa,10);
						throw new RuntimeException("previous atom position was still occupied");
						
					}
					*/ 
					enforce_atom(enh_pix,impWidth,positions[a],weights[a]);
					if( single_atom )
					{
						//condensate if a single atom is nudged towards another one
						atomSolidD = atomSolidD_backup;
						put_down_atom(enh_pix,impWidth,positions[a],Double.NaN);
						atomSolidD = atomSolidD_tiny; 
					}
					
					sim_pix = sparse_floor_fold(enh_pix,impWidth, modelFloor,a_pix,atomWidth);
					model_merit = img_match(sim_pix,imp_pix,impWidth);
					push_pixels(0,true);	
					
					
					if(single_atom && (!protect_atoms))
					{	continue;} //skip weight adjustments altogether
					
					//step2 tune the intensities
					if(optimize_weights)
					{	
						posa = positions[a];
						weighta = weights[a];
						if( (posa == -1) || (enh_pix[posa] == 0.0f) ) //some other atom has eaten this one
						{
							positions[a] = -1;
							continue;
						}
						/*if(false) //quick and dirty
						{
							weights[a] = fit_local_weight(enh_pix, posa, weighta);
						}else*/
						{
							double f = 0.8; //hmm how far off may weights be?
							double fmax = 0.99;
							boolean shrinking = random.nextBoolean();
							while( (weights[a] != 0.0) && f < fmax)
							{
								final double t_weight = shrinking ? f*weights[a] : weights[a]/f;
								enforce_atom(enh_pix,impWidth,posa,t_weight);
								sim_pix = sparse_floor_fold(enh_pix,impWidth, modelFloor,a_pix,atomWidth);
								final double t_merit = img_match(sim_pix,imp_pix,impWidth);
								if(t_merit > model_merit)
								{
									model_merit = t_merit;
									weights[a] = (float)t_weight;
								}
								else
								{
									shrinking = !shrinking;
									f=Math.sqrt(f);
									//next or final put down will cancel the wrong weight anyways
								}	
							}
							
						}
						//kill optimized and weak even when atoms are protected
						if(  (weights[a] < peakWeightMin * avg_peak_weight) && (weights[a] > peakWeightMin * avg_hole_weight))
						{	weights[a] = 0.0f;}
						enforce_atom(enh_pix,impWidth,positions[a],weights[a]);//takes also care of modelFloor
						sim_pix = sparse_floor_fold(enh_pix,impWidth, modelFloor,a_pix,atomWidth);
						push_pixels(0,true);
						model_merit = img_match(sim_pix,imp_pix,impWidth);
					}
					
					if(!single_atom)
					{
						push_pixels(0,true);
						if(seat == 0)
						{	IJ.showProgress( fr*nA+a, full_runs*nA);}
					}	
					
					
					if( (!single_atom) && IJ.escapePressed())
					{
						endless = false;
						cancelled = true;
					}
					
				}
				
				
					
				//step 3 eliminate collided atoms
				//skip this step int last runs if atoms are protected, leave the choice over collisions to the user
				//also skip it if there was no nudging anyways
				final int rmax = (int)Math.ceil(wobble);
				if(!single_atom && (rmax > 0) && ((fr < full_runs-1) || !protect_atoms) )
				{
					atomSolidD = atomSolidD_backup; //re-extend deathzones "activate traps"
					//System.out.println("Before condesation: Atoms: " + numAtoms + " Holes: " + numHoles + " modelFloor: " + modelFloor +
					//						" avg_atom: " + avg_peak_weight + " avg_hole: " + avg_hole_weight );						
					//						+ " impTotal: " + impTotal + " validArea: " + validArea
					/*
					System.out.println("calculated: " + ( ((impTotal-(total_peak_weight+total_hole_weight))/validArea )) );						
					*/
					//check_floor(true);
					for(int i = 1; i < nA; ++i)
					{
						final int j = random.nextInt(i+1);
						if(i!=j)
						{
							final int p = positions[j];
							positions[j] = positions[i];
							positions[i] = p;
							final float w = weights[j];
							weights[j] = weights[i];
							weights[i] = w;		
						}
					}
					
					for(int a = 0; a < nA; ++a)
					{	
						int posa = positions[a];
						float weighta = weights[a];
						//check ifsome other atom has crashed into this one even with tiny solid Diameter
						if( (posa == -1) || (enh_pix[posa] == 0.0f) ) 
						{
							continue;
						}
						/*
						if(enh_pix[posa]!=weights[a]) //passes if we dont use NaN instead of weight
						{
							throw new RuntimeException("weight[a] and enh_pix[posa] dont match up a: " + a);
						}
						*/ 
						//NaN signals to clear all contributions within solid disk
						//and condenses their sum at the biggest atom/hole
						put_down_atom(enh_pix,impWidth,posa,Double.NaN);
					}
					if(fr < full_runs-1)
					{
						//IJ.log("Condensation step, atoms: " + numAtoms + " holes: " + numHoles);
						init_mesh(false);
						show_mesh(mesh,enhanced);
					}
					
					//System.out.println("After condensation: Atoms: " + numAtoms + " Holes: " + numHoles + " modelFloor: " + modelFloor +
					//						" avg_atom: " + avg_peak_weight + " avg_hole: " + avg_hole_weight);
					// 						+ " impTotal: " + impTotal + " validArea: " + validArea
					/*
					System.out.println("calculated: " + ( ((impTotal-(total_peak_weight+total_hole_weight))/validArea )) );
					*/
				}	
					
					
				
				if( (!single_atom) && IJ.escapePressed())
				{
					endless = false;
					cancelled = true;
				}
				sim_pix = sparse_floor_fold(enh_pix,impWidth, modelFloor,a_pix,atomWidth);	
				push_pixels(0,true);
				//check_floor(true);
			} // end for runs
			if( (!single_atom) && (seat == 0) && (!do_stack) )
			{	
				IJ.log("thread: " + seat + " finished " + full_runs + " tuning on " +
				 impT + "(" +impSlice+ ") merit: " + model_merit + " atoms: " + numAtoms + " holes: " + numHoles);
			}
		}
		while(endless && (!cancelled) );
		atomSolidD = atomSolidD_backup;//just in case we skipped final condensation
		return positions[0];
	}
	
	void pull_pixels(int z)
	{
		if(z==0)
		{z=impSlice;}
		
		a_pix = (float[])atomSt.getPixels(1);
		imp_pix = (float[])impSt.getPixels(z);
		if(inverted)
		{
			imp_pix = Arrays.copyOf(imp_pix,imp_pix.length);
			for(int i = 0; i < imp_pix.length; ++i)
			{	imp_pix[i] = -imp_pix[i];}
		}
		
		sim_pix = (float[])simulatedSt.getPixels(z);
		enh_pix = (float[])enhancedSt.getPixels(z);
		dev_pix = (float[])deviationsSt.getPixels(z);
		
	}
	
	void push_pixels(int z,boolean force_dev)
	{
		if(z==0)
		{z=impSlice;}
		
		atomSt.setPixels(a_pix,1);
		if( (atom==null) || (!atom.isVisible()) && (nextHelperSlice.get()!=0))
		{
			nextHelperSlice.set(0);
			dialogs.clear();
			IJ.log("thread: " + seat + " canceled all remaining runs & tasks");
			//System.out.println("thread:" + seat + " cancelled all remaining runs");
		}
		
		simulatedSt.setPixels(sim_pix,z);
		if(simulated.getSlice()==z)
		{	
			//simulated.setStack(simulatedSt);
			simulated.updateAndRepaintWindow();
		}
		
		enhancedSt.setPixels(enh_pix,z);
		//enhanced.setStack(enhancedSt);
		if(enhanced.getSlice()==z)
		{	
			//enhanced.setStack(enhancedSt);
			enhanced.updateAndRepaintWindow();
		}
		
		if(force_dev || deviations.getSlice()==z)
		{
			if(!hex_pixels)
			{
				for(int ind = 0; ind < impArea; ++ind)
				{	
					final float ival = imp_pix[ind];
					//dev_pix[ind] = Float.isNaN(ival)?0.0:ival-sim_pix[ind];
					dev_pix[ind] = ival-sim_pix[ind];
				}
			}
			else
			{
				for(int h = 0; h < imgHP.hpMax; ++h)
				{	
					final int ind = imgHP.hp[h];
					dev_pix[ind] = imp_pix[ind]-sim_pix[ind];
				}
			}
			deviationsSt.setPixels(dev_pix,z);
			if(deviations.getSlice()==z)
			{	
				//deviations.setStack(deviationsSt);
				deviations.updateAndRepaintWindow();
			}
		}
	}
	
	public boolean doDialog()
	{
		if( dialogs.size() == 0) //no pending tasks, get more form user
		{
			do
			{	
				NonBlockingGenericDialog gd = null;
				do
				{
					gd = new NonBlockingGenericDialog(getClass().getSimpleName() + 
								((external_merit_source!=null)?" with "+external_merit_source:""));
					if(!getInput(gd))
					{	return false;}//user cancelled all tasks 
				}
				while(readDialog(gd,true) == true && validateInput(false) != 0); //false = init_stuff
				//System.out.println("sheduling task: " + dialogs.size());
				dialogs.addLast(gd);	
			}
			while(append_another_task);	
		}
		//System.out.println("commencing next task, remaining: " + (dialogs.size()-1));
		return readDialog(dialogs.pollFirst(),false); // no peeking actually consume the reader and disable a few options
	}
	
	public boolean getInput(NonBlockingGenericDialog gd)
	{
		
		reset_enhanced = false;//really no need to repeat that by default
		int[] idArray = WindowManager.getIDList();
		int idlen = 0;
		if(idArray != null)
		{	idlen = idArray.length;}
		String[] imptitleArray = new String[idlen+1];
		int extra_atoms = WindowManager.getImage(atomT)==null?1:0;
		String[] atomtitleArray = new String[idlen+extra_atoms];
		
		int extra_sims = WindowManager.getImage(simulatedT)==null?1:0;
		String[] simulatedtitleArray = new String[idlen+extra_sims];
		
		int extra_model = WindowManager.getImage("model.tif")==null?1:0;
		int extra_dual_model = WindowManager.getImage("dual_model.tif")==null?1:0;
		String[] enhancedtitleArray = new String[idlen+extra_model + extra_dual_model];
		
		int extra_devs = WindowManager.getImage(deviationsT)==null?1:0;
		String[] deviationstitleArray = new String[idlen+extra_devs];
		
		if(extra_atoms == 1)
		{	atomtitleArray[0] = atomT;}
		if(extra_sims == 1)
		{	simulatedtitleArray[0] = simulatedT;} 
		
		{
			int n = 0;
			if(extra_model == 1)
			{	enhancedtitleArray[n++] = "model.tif";}
			if(extra_dual_model == 1)
			{	enhancedtitleArray[n++] = "dual_model.tif";} 
		}
		if(extra_devs == 1)
		{	deviationstitleArray[0] = deviationsT;}
		imptitleArray[idlen] = "<current>";
		for (int i = 0; i < idlen; ++i)
		{	
			String title = WindowManager.getImage(idArray[i]).getTitle();
			imptitleArray[i] = title;
			atomtitleArray[i+extra_atoms] = title; 
			simulatedtitleArray[i+extra_sims] = title;
			enhancedtitleArray[i+extra_model+extra_dual_model] = title;
			deviationstitleArray[i+extra_devs] = title;
		}
		
		gd.setSmartRecording(true);
		gd.addMessage("Task: " + dialogs.size());
		
		gd.addChoice("raw image (GRAY32)", imptitleArray, impT);
		
		gd.addCheckbox("inverted",inverted);
		gd.addCheckbox("hex pixels & periodic domains", hex_pixels);
		gd.addCheckbox("periodic image", periodic);
		
		gd.addCheckbox("full stack", do_stack);
		gd.addNumericField("cycles (<0 forever)", cycles, 0);
		gd.addNumericField("max threads", maxThreads, 0);
		
		gd.addMessage("The atom will be reused and overrule the here set FWHM");
		
		gd.addChoice("atom (psf)", atomtitleArray, atomT);
		gd.addNumericField("solid diameter", atomSolidD0, 2);
		gd.addNumericField("guessed atom FWHM", atomFWHM, 2);
		
		gd.addCheckbox("interactive atom", do_setup_atom);
		gd.addCheckbox("enable holes", enable_holes);
		gd.addCheckbox("enforce solid diameter",enforce_solid);
		gd.addMessage("all other images will be reset if dimensions dont match");
		
		gd.addChoice("model", enhancedtitleArray, enhancedT);
		gd.addCheckbox("reset model", reset_enhanced);
		gd.addCheckbox("homogenize intensities", equalize_intensities);
		
		gd.addNumericField("intensity mixing factor, 0 .. enforce treshold",mix,2);
		gd.addCheckbox("fixed intensities", !optimize_weights);
		gd.addCheckbox("generate atoms", optimize_enhanced);
		gd.addCheckbox("and keep all atoms", protect_atoms);
		
		gd.addNumericField("relative peak weight treshold(0..1)", peakWeightMin,2);
		gd.addNumericField("grand trials (<0 one atom)", grand_runs, 0);
		gd.addCheckbox("analyze & fix topology",fix_topology);
		
		gd.addCheckbox("show topology overlay", show_topology);
		gd.addCheckbox("fine tune atoms", optimize_fine);
		gd.addNumericField("wobble radius (>=0)", wobble, 2);
		
		gd.addNumericField("cycles per atom(<0 forever)", fine_runs, 0);
		gd.addCheckbox("optimize_FWHM", optimize_FWHM);
		gd.addCheckbox("auto scaled intensities", optimize_Brightness);
		
		gd.addCheckbox("show dual mesh",count_atoms);
		gd.addMessage("These will be always updated");
		gd.addChoice("simulated", simulatedtitleArray, simulatedT);
		
		gd.addChoice("deviations", deviationstitleArray, deviationsT);
		gd.addCheckbox("init mesh from xyz", mesh_from_xyz);
		gd.addCheckbox("read signature from disk", read_sig_file);
		gd.addCheckbox("join meshes (after all cyles)",join_meshes);
		
		gd.addCheckbox("write signature to disk", write_sig_file);
		gd.addCheckbox("write topology to disk", write_top_file);
		gd.addCheckbox("append another task", append_another_task);
				
		if(loc != null)
		{
			gd.centerDialog(false);
			gd.setLocation(loc.x+10,loc.y+10); //compensate the drift
		}
		if(rec != null)
		{
			gd.setSize(rec.width,rec.height);
		}
		gd.showDialog();
		loc = gd.getLocation();
		rec = gd.getBounds();
		if(gd.wasCanceled())
		{
			return false;
		}
		return true;
	}
	
	public boolean readDialog(NonBlockingGenericDialog gd, boolean peek)
	{	
			try
			{
			Vector<Choice> choices = (Vector<Choice>)gd.getChoices();
			Vector<TextField> numbers = (Vector<TextField>)gd.getNumericFields();
			Vector<Checkbox> boxes = (Vector<Checkbox>)gd.getCheckboxes();
			/*
			System.out.println(choices.toString());
			System.out.println("");
			System.out.println(numbers.toString());
			System.out.println("");
			System.out.println(boxes.toString());
			*/
			int c = 0; //counter for Strings from Choices
			int n = 0; //counter for numbers from TextFields inputs
			int b = 0; //counter for booleans from Checkboxes
			
			impT = choices.get(c++).getSelectedItem(); //there would also be getSelectedIndex()
			
			inverted = boxes.get(b++).getState(); 
			hex_pixels = boxes.get(b++).getState(); 
			periodic = boxes.get(b++).getState();
			
			do_stack = boxes.get(b++).getState();
			cycles = (int)Double.parseDouble(numbers.get(n++).getText());
			maxThreads = (int)Double.parseDouble(numbers.get(n++).getText());
		
			atomT = choices.get(c++).getSelectedItem();
			atomSolidD0 = Double.parseDouble(numbers.get(n++).getText());
			atomFWHM = Double.parseDouble(numbers.get(n++).getText());
		
			do_setup_atom = boxes.get(b++).getState() && peek;//only effective upon peeking
			enable_holes = boxes.get(b++).getState();
			enforce_solid = boxes.get(b++).getState();
			enhancedT = choices.get(c++).getSelectedItem();
			
			reset_enhanced = boxes.get(b++).getState();
			equalize_intensities = boxes.get(b++).getState();
			mix = Double.parseDouble(numbers.get(n++).getText());
			
			optimize_weights = !boxes.get(b++).getState();
			optimize_enhanced = boxes.get(b++).getState();
			protect_atoms = boxes.get(b++).getState();
			peakWeightMin = Double.parseDouble(numbers.get(n++).getText());
			
			grand_runs = (int)Double.parseDouble(numbers.get(n++).getText());
			fix_topology = boxes.get(b++).getState();
			show_topology = boxes.get(b++).getState();
		
			optimize_fine = boxes.get(b++).getState();
			wobble = Math.abs(Double.parseDouble(numbers.get(n++).getText()));
			fine_runs = (int)Double.parseDouble(numbers.get(n++).getText());
			
			optimize_FWHM = boxes.get(b++).getState();
			optimize_Brightness = boxes.get(b++).getState();
			count_atoms = boxes.get(b++).getState();
			
			simulatedT = choices.get(c++).getSelectedItem();
			deviationsT = choices.get(c++).getSelectedItem();
			mesh_from_xyz = boxes.get(b++).getState();
			read_sig_file = boxes.get(b++).getState();
			
			join_meshes = boxes.get(b++).getState();
			write_sig_file = boxes.get(b++).getState();
			write_top_file = boxes.get(b++).getState();
			
			append_another_task = boxes.get(b++).getState() && peek; //only effective upon peeking
			}
			catch(Exception e)
			{
				e.printStackTrace();
				return false;
			}
		/*
		//keep this block for debugging
		{
			impT = gd.getNextChoice();
			hex_pixels = gd.getNextBoolean();
			periodic = gd.getNextBoolean();
			
			do_stack = gd.getNextBoolean();
			cycles = (int)gd.getNextNumber();
			maxThreads = (int)gd.getNextNumber();
			
			atomT = gd.getNextChoice();
			atomSolidD0 = gd.getNextNumber();
			atomFWHM = gd.getNextNumber();
			
			do_setup_atom = gd.getNextBoolean();
			do_setup_atom = false; //this option has already bee processed by "dry" validation
			enforce_solid = gd.getNextBoolean();
			enhancedT = gd.getNextChoice();
			
			reset_enhanced = gd.getNextBoolean();
			equalize_intensities = gd.getNextBoolean();
			mix = gd.getNextNumber();
			
			optimize_weights = !gd.getNextBoolean();
			optimize_enhanced = gd.getNextBoolean();
			peakWeightMin = gd.getNextNumber();
		
			grand_runs = (int)gd.getNextNumber();
			fix_topology = gd.getNextBoolean();
			show_topology = gd.getNextBoolean();
			
			optimize_fine = gd.getNextBoolean();
			wobble = Math.abs(gd.getNextNumber());
			fine_runs = (int)gd.getNextNumber();
			
			optimize_FWHM = gd.getNextBoolean();
			optimize_Brightness = gd.getNextBoolean();
			count_atoms = gd.getNextBoolean();
			
			simulatedT = gd.getNextChoice();
			deviationsT = gd.getNextChoice();
			mesh_from_xyz = gd.getNextBoolean();
			read_sig_file = gd.getNextBoolean();
			
			join_meshes = gd.getNextBoolean();
			write_sig_file = gd.getNextBoolean();
			write_top_file = gd.getNextBoolean();
			
			append_another_task = gd.getNextBoolean();
			append_another_task = false; //this is just not allowed when not peeking
		}
		*/
		return true;
	}
	
	
	int validateInput(boolean init_stuff) //0..ok, 1 .. auto fixed input, 2 .. broken input
	{
		if(!init_stuff)
		{
			IJ.log("sheduling Task: " + dialogs.size());
		}
		int ic = instance_count.get();
		int res = 0;
		if(hex_pixels)
		{
			if(periodic)
			{
				IJ.log("hexagonal pixels are always periodic inside hex domains, disabling torus periodicity");
				periodic = false;
				return 2;
			}
			if(read_sig_file || join_meshes || write_sig_file || write_top_file || mesh_from_xyz)
			{
				read_sig_file = false;
				join_meshes = false;
				write_sig_file = false;
				write_top_file = false;
				mesh_from_xyz = false;
				IJ.log("hex_pixels cannot be combined with anything that requires Signatures");	
				return 2;
			}
		}
		
		if(ic > 1)
		{
			IJ.showMessage("There are still "+ (ic-1) +" threads running, try again later");
			return 2;
		}
		if(enhancedT.equals("model.tif") && inverted )
		{
			IJ.log("model.tif is a reserved default name for direct lattices and regular images");
			return 2;
		}
		if(enhancedT.equals("dual_model.tif") && (!inverted) )
		{
			IJ.log("dual_model.tif is a reserved default name for dual lattices and inverted images");
			return 2;
		}
		
		if(maxThreads < 1)
		{	
			IJ.log("Cannot employ at most " + maxThreads);
			return 2;	
		}
		if(external_merit) //merit is supposed to come from reconstruction
		{
			if(maxThreads != 1)
			{
				IJ.log("There can only be one single thread if an external merit is to be used");
				maxThreads = 1;
				return 2;
			}
			if(!hex_pixels)
			{
				IJ.log("The external merit from " + external_merit_source + " does require hexagonal pixels");
				hex_pixels = true;
				return 2;
			}
		}
		
		
		if( ( ( ((!init_stuff) || do_stack) && (grand_runs < 0) || (fine_runs < 0) ) ) ||
			( (grand_runs < 0) && (fine_runs < 0) ) 	)
		{
			IJ.log("impossible request?");
			IJ.log("another task: " + (!init_stuff));
			IJ.log("cycles: " + cycles);
			IJ.log("do stack: " + do_stack);
			IJ.log("grand runs: " + grand_runs);
			IJ.log("fine_runs: " + fine_runs);
			if(!IJ.showMessageWithCancel("Seriously?","Proceed at your own risk!"))
			{	return 2;}
		}
		if(impT.equals("<current>"))
		{	
			imp = IJ.getImage();
			if(imp != null)
			{	impT = imp.getTitle();}	
		}
		else
		{	imp = WindowManager.getImage(impT);}
		atom = WindowManager.getImage(atomT);
		enhanced = WindowManager.getImage(enhancedT);
		simulated = WindowManager.getImage(simulatedT);
		deviations = WindowManager.getImage(deviationsT);
		if(imp == null)
		{
			IJ.log("Error: could not find " + impT);
			return 2; 
		}
		else if(imp.getType() != ImagePlus.GRAY32 )
		{
			IJ.log("Error: only float images are supported, please convert");
			return 2; //not yet a real issue
		}
		if( atom != null )
		{
			atomWidth = atom.getWidth();
			atomHeight = atom.getHeight();
			if(hex_pixels)
			{
				if((atomWidth != atomHeight) || (atomWidth%2 == 1) )
				{
					IJ.log("hex pixels do require square images with even edge lengths");
					do_setup_atom = true;
					return 2;
				}
				atomHP = new HexPixels(atomWidth/2);
			}
			else
			{	atomHP = null;}
			measure_FWHM((float[])atom.getStack().getPixels(1));
		
		}
		
		
		if( ( (atom == null) || do_setup_atom) && (!setup_atom()) )
		{	
			IJ.log("nothing to do without " + atomT);
			do_setup_atom = true;
			return 2;	
		}
		
		atomWidth = atom.getWidth();
		atomHeight = atom.getHeight();
		if( (!hex_pixels) && ((atomWidth % 2 == 0) || (atomHeight % 2 == 0)) )
		{
			IJ.log("Error: for square pixels " + atomT +" must be a kernel with odd edge lengths");
			do_setup_atom = true;
			return 2;
		}
		if( (hex_pixels) && ((atomWidth % 2 == 1) || (atomWidth != atomHeight )) )
		{
			IJ.log("Error: for hex pixels and domain " + atomT +" must be a square with even edge length");
			do_setup_atom = true;
			return 2;
		}
		
		if(atomSolidD0 < 2)
		{
			IJ.log("solid diameter of atom cannot be " + atomSolidD0);
			do_setup_atom = true;
			return 2;
		}
		do_setup_atom = false;//dont repeat that
		atomSolidD = atomSolidD0;
		impWidth = imp.getWidth();
		impHeight = imp.getHeight();
		impDepth = imp.getStackSize();
		impSlice = imp.getSlice();
		impSlice0 = impSlice;
		impArea = impWidth * impHeight;
		if(hex_pixels)
		{
			if( (impWidth != impHeight) || (impWidth % 2 == 1))
			{
				IJ.log("Error: for hex pixels and domain " + impT +" must be a square with even edge length");
				return 2;
			}
			if(init_stuff)
			{
				if( (atomHP == null) || (atomHP.get_Width() != atomWidth))
				{	atomHP = new HexPixels(atomWidth/2);}
				imgHP = new HexPixels(impWidth/2);
				tinyHP = new HexPixels((int)(atomSolidD+0.5));
				boxHP = new HexPixels(atomWidth); //special case for final int r = 1;
			}
		}
		else
		{
			atomHP = null;
			imgHP = null;
			tinyHP = null;
			boxHP = null;
		}
		
		if(init_stuff)
		{
			//set the defaults so that SurfaceMesh_Options can actually show them
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
			Smart_Atomizer.master_mesh = new SurfaceMesh(impWidth,impHeight,0);
			SurfaceMesh.master_mesh = Smart_Atomizer.master_mesh;
			master_mesh.hex_pixels = hex_pixels;
			master_mesh.imgHP = imgHP;
			master_mesh.periodic = periodic;
			master_mesh.inverted = inverted;
			master_sig_string = read_sig_file ? read_sig() : null;
			if(master_sig_string != null)
			{	master_mesh.link_signature(master_sig_string);}
		}
		
		numTasks = (do_stack&&(!optimize_FWHM) ? impDepth : 1); //we can have at most one thread per slice
		if(!(optimize_enhanced || optimize_fine)) //no multithreading
		{ 
			//numTasks = 1;
			if(cycles < 0)
			{
				IJ.log("refusing to repeat trivial task forever");
				cycles = 1;
			}
			
		}
		if((!init_stuff) && do_setup_atom)
		{	
			IJ.log("There cannot be interactive atoms for prescheduled tasks");
			do_setup_atom = false;
			return 2;
		}
		
		
		if(numTasks > taken_slices.length())//but we have to respect the hard coded limit
		{	numTasks = taken_slices.length();}
		if( (maxThreads > 0) && (numTasks > maxThreads) )
		{	numTasks = maxThreads;} //limit threads to user set number
		/* //No longer an issue as the last remaining thread will become the new master
		if( (numTasks > 1) && (cycles > 1) && (impDepth <= numTasks) )
		{	numTasks = impDepth-1;} //dont run all slices at once if a certain number of cylcles is requested
		*/
		final int total_runs = (do_stack ? impDepth : 1) * cycles;
		if(do_stack && optimize_FWHM)
		{	
			IJ.log("Cannot allow multi threading for adjusting FWHM");
			return 2;
		}
		
		if(numTasks > 1 && mesh_from_xyz)
		{
			IJ.log("multithreading and manually choosing xyz files is disabled");
			return 2;
		}
		
		
		if(!init_stuff)
		{
			if(do_setup_atom)
			{	
				IJ.log("There cannot be interactive atoms for prescheduled tasks");
				do_setup_atom = false;
				return 2;
			}
			IJ.log("inv. image: " + inverted + " holes: " + enable_holes + " hex_pixels: " + hex_pixels + " periodic: " + periodic + " full stack: " + do_stack + " cycles: " + cycles + " total runs: " + total_runs + " max threads: " + numTasks);
			IJ.log("relative peak weight threshhold: " + peakWeightMin + " atom solid diameter: " + atomSolidD0 + " preserving atoms: " + protect_atoms);
			IJ.log("equalize intensities: " +  equalize_intensities + " with mixing: " + mix + " auto Brightness: " + optimize_Brightness);
			IJ.log((optimize_enhanced?("grand_runs: " + grand_runs):"") + " fix_mesh: " + fix_topology + (optimize_fine?" fine_runs: " + fine_runs:"") +
			" wobble: " + wobble + " const_weights: " + (!optimize_weights) );
			
			
			
			IJ.log("**************Settings for Task: " + dialogs.size() + " passed validation*****************");
		}
		/*else
		{
			IJ.log("Vou may rerun " + getClass().getSimpleName() + " to cancel pending cycles" );
		}*/
		if(init_stuff)
		{
			for(int d=0; d < taken_slices.length(); ++d)
			{
				taken_slices.set(d,0);//this thread is supposed to be idle
			}
			ensure_all_images(impSlice);
			very_first_slice = true; //the master will not have to repeat this
			impSt = imp.getImageStack();
			atomSt = atom.getImageStack();
			enhancedSt = enhanced.getImageStack();
			simulatedSt = simulated.getImageStack();
			deviationsSt = deviations.getImageStack(); 
		}
		return res;	
	}
	
	void ensure_all_images(int nextSlice)
	{
		
		impSlice = nextSlice;
		if(seat == 0)
		{	imp.setSlice(nextSlice);}
		
		validArea = impArea;
		{
			float[] pix = (float[])imp.getStack().getPixels(impSlice);
			impAvg = get_avg(pix,imgHP); //also sets global validArea !!!
			impStd = Math.sqrt(get_var(pix,(float)impAvg,imgHP));
			if(inverted)
			{	impAvg = -impAvg;}
			//IJ.log(impT+"("+ impSlice  + ") avg: " + impAvg + ", std: " + impStd);
		}
		mesh = meshes[impSlice-1];
		impTotal = impAvg*(hex_pixels ? imgHP.hpMax : validArea);
		enhanced = ensure_image(enhanced,enhancedT,impWidth,impHeight,impDepth);
		simulated = ensure_image(simulated,simulatedT,impWidth,impHeight,impDepth);
		deviations = ensure_image(deviations,deviationsT,impWidth,impHeight,impDepth);
	}
	
	ImagePlus ensure_image(ImagePlus nimp, String nimpT, int width, int height, int depth)
	{
		
		if( reset_enhanced && (nimp == enhanced) && (nimp != null) &&
			(nimp.getWidth() == width) && 
			(nimp.getHeight() == height) &&
			(nimp.getStackSize() == depth))
		{
			ImageStack nimSt = nimp.getStack();
			float[] npix = (float[])nimSt.getPixels(impSlice);
			init_model(npix);
			if(seat==0)
			{
				nimp.setSlice(impSlice);
				nimp.show();
			}
			return nimp;
		}
		
		if(	(nimp == null) || (nimp.getWidth() != width) || 
			(nimp.getHeight() != height) || (nimp.getStackSize() != depth) )
		{
			if(nimp != null)
			{	
				IJ.log("dismissing outdated " + nimpT);
				nimp.close();
			}
			nimp = NewImage.createFloatImage(nimpT,width,height,depth,NewImage.FILL_BLACK);
			//IJ.log("created " + nimpT);
			if(nimpT.equals(enhancedT))
			{
				nimp = init_enhanced(nimp, 0 );
			}
			
		}
		if(nimpT.equals(simulatedT))
		{
			nimp = init_simulated(nimp, 0 );
		}
		if(nimpT.equals(deviationsT))
		{
			nimp = init_deviations(nimp, 0 );
		}	
		if(seat == 0)
		{
			nimp.setSlice(impSlice);
			nimp.show();
		}
		return nimp;
	}
	
	ImagePlus init_enhanced(ImagePlus nimp, int z)
	{
		if(enhanced != null)
		{
			IJ.log("Error, cannot initialize " + enhancedT + " with existing " + enhancedT);
			return nimp;
		}
		if(imp == null)
		{
			IJ.log("Error, cannot make initial guess for " + enhancedT + " without " + impT);
			return nimp;
		}
		
		//ImageStack imSt = imp.getStack();
		if(z==0) //requested slice
		{ z = impSlice;}
		if(seat == 0)
		{	IJ.showStatus("updating " + enhancedT + "("+ z + ")");}
		float[] npix = new float[impArea];
		
		init_model(npix);
		
		ImageStack nimpSt = nimp.getStack();
		nimpSt.setPixels(npix,z);
		nimp.setStack(nimpSt);
		return nimp;
	}
	
	ImagePlus init_simulated(ImagePlus nimp, int z )
	{
		if(enhanced == null)
		{
			IJ.log("Error, cannot initialize " + simulatedT + " without " + enhancedT);
			return nimp;
		}
		if(atom == null)
		{
			IJ.log("Error, cannot initialize " + simulatedT + " without " + atomT);
			return nimp;
		}
		
		ImageStack enSt = enhanced.getStack();
		if(z==0) //requested slice
		{ z = impSlice;}
		if(seat == 0)
		{	IJ.showStatus("updating " + simulatedT + "("+ z + ")");}
		float[] enpix = (float[])enSt.getPixels(z);
		final int area = enpix.length;
		if(enforce_solid)
		{	enforce_solid_diameter(enpix);}
		
		count_and_report_atoms(enpix,(seat==0));
		
		
		if( (numAtoms + numHoles > 0) && equalize_intensities)
		{
			if(mix >= 0 && mix <= 1.0)
			{
				final float peak_lvl = (float) (peakWeightMin * avg_peak_weight);
				final float hole_lvl = (float) (peakWeightMin * avg_hole_weight);
				if(seat == 0)
				{
					if(mix > 0.0)
					{	IJ.log("equalizing intensities in " + enhancedT + " with mixing: " + mix);}
					else
					{	IJ.log("tresholding intensities at " + peakWeightMin + " * mean : " + peak_lvl + "/" + hole_lvl);}
				}
				for(int ind = 0; ind < area; ++ind)
				{
					final float enval = enpix[ind];
					if(enval > 0.0f)
					{
						if(mix > 0.0)
						{	
							enpix[ind] = (float)((1.0-mix)*enval + mix*avg_peak_weight);
						}
						else if (mix == 0.0)
						{
							if(enval < peak_lvl )
							enpix[ind] = 0.0f;
						}
					}
					else if(enval < 0.0f)
					{
						if(mix > 0.0)
						{	
							enpix[ind] = (float)((1.0-mix)*enval + mix*avg_hole_weight);
						}
						else if (mix == 0.0)
						{
							if(enval > hole_lvl )
							enpix[ind] = 0.0f;
						}
					}
				}
			}
			if(mix==0.0f)
			{	count_and_report_atoms( enpix, (seat==0) );}
		}
		ImageStack atSt = atom.getStack();
		float[] atompix = (float[])atSt.getPixels(1);
		float[] npix = sparse_floor_fold( enpix, impWidth, modelFloor,  atompix, atomWidth);
		ImageStack nimpSt = nimp.getStack();
		nimpSt.setPixels(npix,z);
		nimp.setStack(nimpSt);
		return nimp;
	}
	
	void count_and_report_atoms(float[] enpix, boolean report)
	{	
		final int area = hex_pixels?imgHP.hpMax:enpix.length;
		double ensum = 0.0;
		double ensump = 0.0;
		double ensumm = 0.0;
		int countp = 0;
		int countm = 0;
		for(int ind = 0; ind < area; ++ind)
		{
			final int i = hex_pixels?imgHP.hp[ind]:ind;
			final float val = enpix[i];
			if( (val > 0.0f) || (val < 0.0f) )
			{
				if(val > 0.0f)
				{	
					++countp;
					ensump += val;
				}
				else if(val < 0.0f)
				{	
					++countm;
					ensumm += val;	
				}
				
				ensum += val;
				
			}
		}
		total_peak_weight = ensump;
		total_hole_weight = ensumm;
		numAtoms = countp;
		numHoles = countm;
		final double lastfloor = modelFloor;
		//if(inverted)
		//{	ensum = -ensum;}
		modelFloor = ((impTotal-ensum)/(hex_pixels ? imgHP.hpMax : validArea));
		avg_peak_weight = ( (numAtoms > 0) ? total_peak_weight/numAtoms : 0.0 );
		avg_hole_weight = ( (numHoles > 0) ? total_hole_weight/numHoles : 0.0 );
		final int coverage = (int)(100*(total_peak_weight-total_hole_weight)/impTotal);
		if(report)
		{	IJ.log("thread: "+seat+" finds " + numAtoms + " atoms and "+ numHoles + " holes in "+ enhancedT+ "("+ impSlice +
			") accounting for " + coverage + "% of the total intensity. avg. int.: " + avg_peak_weight + "/" + avg_hole_weight);
			//IJ.log("last modelFloor: "+ lastfloor +" new modelFloor: " + modelFloor);
		}
		//System.out.println("ensum: " + ensum + " impTotal: " + impTotal + " simsum: " + (ensum+area*modelFloor) + " modelFloor: " + modelFloor);
		//check_floor(true); //the test would be run on the previous model!
	}
	
	void enforce_solid_diameter(float[] enpix)
	{
		float[] ipix = (float[])imp.getStack().getPixels(impSlice);
		int filtered = 0;
		for(int ind = 0; ind < (hex_pixels ? imgHP.hp.length : enpix.length) ; ++ind)
		{
			final int i = hex_pixels?imgHP.hp[ind]:ind;
			final float val = enpix[i];
			if( (val > 0.0f) || (val < 0.0f) )
			{
				/*
				if(protect_atoms)
				{
					enforce_atom(enpix, impWidth, i, Float.isNaN(ipix[i])?0.0f:val);
				}
				else
				{
					
				}
				*/ 
				put_down_atom(enpix, impWidth, i, Float.isNaN(ipix[i])?0.0f:Double.NaN); //trigger condensation
				filtered += deltaAtoms;
			}
		}
		if(seat==0)
		{	IJ.log("enforced solid diameter eliminated atoms: " + filtered);}
	}
	
	
	ImagePlus init_deviations(ImagePlus nimp, int z )
	{
		if(imp == null)
		{
			IJ.log("Error, cannot calculate " + deviationsT + " without " + impT);
			return nimp;
		}
		if(simulated == null)
		{
			IJ.log("Error, cannot calculate " + deviationsT + " without " + simulatedT);
			return nimp;
		}
		if(z==0) //requested slice
		{ z = impSlice;}
		if(seat == 0)
		{	IJ.showStatus("updating " + deviationsT + "("+ z + ")");}
		ImageStack nimpSt = nimp.getStack();
		float[] npix = (float[])nimpSt.getPixels(z);
		
		ImageStack imSt = imp.getStack();
		float[] impix = (float[])imSt.getPixels(z);
		
		ImageStack simSt = simulated.getStack();
		float[] simpix = (float[])simSt.getPixels(z);
		
		final int area = npix.length;
		for(int ind = 0; ind < area; ++ind)
		{
			npix[ind] = (inverted?-1.0f:1.0f)*impix[ind]-simpix[ind];
		}
		return nimp;
	}
	
	boolean setup_atom()
	{
		int atom_rx = 15;
		int atom_ry = -1;
		boolean accept_atom = false;
		do
		{
			if(atom != null)
			{
				atom_rx = atom.getWidth()/2;
				atom_ry = atom.getHeight()/2;
			}
			
			NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Setup Atom for "+getClass().getSimpleName());
			gd.setSmartRecording(true);
			if(atom != null)
			{
				gd.addCheckbox("accept " + atomT + " as it is", accept_atom);
			}
			gd.addMessage("parameters for next " + atomT);
			gd.addNumericField("radius x",atom_rx,0);
			if(!hex_pixels)
			{	gd.addNumericField("radius y(-1..for ry=rx)",atom_ry,0);}
			
			gd.addNumericField("FWHM", atomFWHM, 2);
			gd.addNumericField("solid diameter", atomSolidD0, 2);
			if(loc != null)
			{
				gd.centerDialog(false);
				gd.setLocation(loc.x+10,loc.y+10); //compensate the drift
			}
			gd.showDialog();
			loc = gd.getLocation();
			if(gd.wasCanceled())
			{
				return false;
			}
			
			accept_atom = gd.getNextBoolean();
			if(accept_atom && (atom != null) ) //second condition is for safety only
			{	return true;}
			atom_rx = (int)gd.getNextNumber();
			if(!hex_pixels)
			{	atom_ry = (int)gd.getNextNumber();}
			else
			{	atom_ry = atom_rx;}
			atomFWHM = gd.getNextNumber();
			atomSolidD0 = gd.getNextNumber();
			
			if(atom_ry == -1)
			{atom_ry = atom_rx;}
			
			if( (atom_rx < 1) || (atom_ry < 1) || (atomFWHM <= 0.0) || atomSolidD <= 0.0)
			{
				IJ.log("stupid user detected, reseting radii and FWHM");
				if(atom_rx < 0)
				{	atom_rx = 15;}
				if(atom_ry < 0)
				{	atom_ry = 15;}
				if(atomFWHM <= 0)
				{	atomFWHM = 10;}
				if(atomSolidD0 <= 0)
				{	atomSolidD0 = atomFWHM;}
				
				continue;
			}
			
			
			atomWidth = 2 * atom_rx + (hex_pixels?0:1);
			atomHeight = 2 * atom_ry + (hex_pixels?0:1);
			atom = ensure_image(atom,atomT,atomWidth,atomHeight,1);
			//inv_atom = ensure_image(inv_atom,inv_atomT,atomWidth,atomHeight,1);
			if(hex_pixels)
			{	atomHP = new HexPixels(atom_rx);}
			init_atom(null);
			
			
		}
		while(true);
	}

	void init_atom(float[] atompix)
	{
		if(atompix == null)
		{
			ImageStack atSt = atom.getStack();
			atompix = (float[])atSt.getPixels(1);
		}
		double sig = atomFWHM / 2.235482;
		double b = 1.0/(2*sig*sig);
		int area = atomWidth * atomHeight;
		
		int px = atomWidth / 2;
		int py = atomHeight / 2;
		final float intensity = 1.0f; 
		if(hex_pixels)
		{
			b*=0.5;
			final int hpMax = atomHP.hpMax;
			final int[] hp = atomHP.hp;
			for(int i = 0; i < hpMax; ++i)
			{
				final int qr = hp[i];
				final int[] r = atomHP.to_xyz(qr);
				float val = (float)Math.exp(-((r[0]*r[0]+r[1]*r[1]+r[2]*r[2])*b));
				atompix[qr] = intensity * val;
			}
		}
		else
		{
		
			for(int ind = 0; ind < area; ++ind )
			{
				int rx = ind % atomWidth - px;
				int ry = ind / atomWidth - py;
				float val = (float)Math.exp(-(rx*rx+ry*ry)*b);
				atompix[ind] = intensity * val;
			}
		}
		normalize(atompix,1.0);
	}

	double normalize(float[] pixels,double norm)
	{
		double sum = 0.0;
		final int area = pixels.length;
		for(int ind = 0; ind < area; ++ind)
		{
			sum += pixels[ind];	
		}
		double asum = Math.abs(sum); //normalize negative images to -1 * norm
		asum = (float)(norm/asum);
		for(int ind = 0; ind < area; ++ind)
		{
			pixels[ind]*=asum;	
		}
		return sum;
	}
	
	//cannot handle hex_pixels, but no longer in use anyways
	float[] foldwith( float[] img, final int imgW, final float[] kernel, final int kerW)
	{
		final int iarea = img.length;
		final int imgH = iarea/imgW;
		final int karea = kernel.length;
		final int kerH = karea/kerW;
		
		final float[] korr = new float[iarea];
		
		for(int i=0;i<iarea;++i)
		{
			final int ix = i%imgW;
			final int iy = i/imgW;
			
			float psum = 0.0f;
			//float pnorm = 0.0f;
			
			for(int k=0;k<karea;++k)
			{
				final int kx = k%kerW;
				final int ky = k/kerW;
				//move origin to center, dont use JAVA % operations for negative indices
				final int kmx = kx - kerW/2;
				final int kmy = ky - kerH/2;	
			
				int px = ix-kmx; //negate here
				int py = iy-kmy; //negate here
				
				if(periodic) //works as long as kernel is smaller than image
				{
					px = (px + imgW)%imgW;
					py = (py + imgH)%imgH;
				}
				
				if( periodic || 
					( (px >= 0) && (px < imgW) &&
					  (py >= 0) && (py < imgH)    )	)
				{
					final float kval = kernel[k];
					final int p = px + py * imgW;
					final float imgval = img[p];
					//pnorm += kval;
					psum += kval * imgval;		
				}
			}
			korr[i] = psum;
		}
		return korr;
	}
	//good implementation for sparse images
	float[] sparse_floor_fold( float[] img, final int imgW, final double dfloor,  final float[] kernel, final int kerW)
	{
		final float mfloor = (float)dfloor;
		final int iarea = img.length;
		final int imgH = iarea/imgW;
		final int karea = kernel.length;
		final int kerH = karea/kerW;
		
		HexPixels imHP =  hex_pixels ? (imgW==impWidth  ? imgHP  : (imgW == 2*atomWidth? boxHP : new HexPixels(imgW/2) )) : null;
		HexPixels kerHP = hex_pixels ? (kerW==atomWidth ? atomHP : new HexPixels(kerW/2) ) : null; 
		
		final float[] korr = new float[iarea];
		
		if(hex_pixels)
		{
			for(int h = 0; h < imHP.hpMax; ++h)
			{	korr[imHP.hp[h]] = mfloor;}
			for(int i = 0; i < imHP.hpMax; ++i)
			{
				final int iqr = imHP.hp[i];
				if(img[iqr]!=0.0f)
				{
					final float ival = img[iqr];
					
					for(int k = 0; k < kerHP.hpMax; ++k)
					{
						final int kqr = kerHP.hp[k];
						final float kval = kernel[kqr];
						final int[] kxyz = kerHP.to_xyz(kqr);
						final int ikqr = imHP.shift_qr(iqr,kxyz);
						korr[ikqr] += ival * kval;
					}
				}
			}
		}
		else
		{
			Arrays.fill(korr,mfloor);
			for(int i=0;i<iarea;++i)
			{
				if(img[i]>0.0f || img[i]<0.0f)
				{
					final int ix = i%imgW;
					final int iy = i/imgW;
					
					//float psum = 0.0f;
					//float pnorm = 0.0f;
					final float ival = img[i];
					
					for(int k=0;k<karea;++k)
					{
						final int kx = k%kerW;
						final int ky = k/kerW;
						//move origin to center, dont use JAVA % operations for negative indices
						final int kmx = kx - kerW/2;
						final int kmy = ky - kerH/2;	
					
						int px = ix+kmx; //positive here
						int py = iy+kmy; //positive here
						
						if(periodic) //works as long as kernel is smaller than image
						{
							px = (px + imgW)%imgW;
							py = (py + imgH)%imgH;
						}
						
						if( periodic || 
							( (px >= 0) && (px < imgW) &&
							  (py >= 0) && (py < imgH)	  ) )
						{
							korr[px+py*imgW] += kernel[k]*ival;
						}	
					}
				}
			}
		}
		if(external_merit)
		{
			for(int i = 0; i < iarea; ++i)
			{
				korr[i] = (float)Math.round(korr[i]);
			}
		}
		return korr;
	}
	
	double get_external_merit(float[] spix)
	{
		try
		{
			while(sharedmerit.get_state() != 0) //Should never be true
			{
				Thread.sleep(10);
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return Double.NaN;
		}
		short[] pixels = (short[])next_modelSt.getPixels(impSlice);
		if(pixels.length != spix.length)
		{
			throw new RuntimeException("length of sim_pix does not match with next_model");
		}
		
		for(int i = 0; i < pixels.length; ++i )	
		{	pixels[i] = (short)sim_pix[i];}
		sharedmerit.set_state(2);//2..requesting update
		return sharedmerit.get_merit();
	}
	
	
	
	double img_match(float[] img1, float[] img2, int imgW)
	{
		double sum1 = 0.0;
		double sum2 = 0.0;
		double std1 = 0.0;
		double std2 = 0.0;
		double mix = 0.0;
		int area = img1.length;
		int varea = 0;
		HexPixels imHP = null;
		if(hex_pixels)
		{
			if(imgW == 2*atomWidth)
			{	imHP = boxHP;}
			else if(imgW == impWidth)
			{	
				imHP = imgHP;
				///FIXME this is to DEBUG the use of external_merit
				 
				if(external_merit)
				{
					// This will and SHOULD crash if the sheredmerit was not provided
					return get_external_merit(img1);
				}
				
			}
			else if(imgW == atomWidth)
			{	imHP = atomHP;}
			else
			{	imHP = new HexPixels(imgW);}
		}
		
		if(!hex_pixels)
		{
			for(int ind = 0; ind < area; ++ind)
			{
				final float v1 = img1[ind];
				final float v2 = img2[ind];
				if( (!Float.isNaN(v1)) && (!Float.isNaN(v2)))
				{
					sum1 += v1;
					sum2 += v2;
					++varea;
				}
			}
		}
		else
		{
			for(int h = 0; h < imHP.hpMax; ++h)
			{
				final int ind = imHP.hp[h];
				sum1 += img1[ind];
				sum2 += img2[ind];
			}
		}
		
		
		final float avg1 = (float)( sum1 / ( hex_pixels ? imHP.hpMax : varea ) );
		final float avg2 = (float)( sum2 / ( hex_pixels ? imHP.hpMax : varea ) );
		
		if(!hex_pixels)
		{
			for(int ind = 0; ind < area; ++ind)
			{
				final float v1 = (img1[ind]-avg1);
				final float v2 = (img2[ind]-avg2);
				if( (!Float.isNaN(v1)) && (!Float.isNaN(v2)) )
				{
					std1 += v1*v1;
					std2 += v2*v2;
					mix += v1*v2; 
				}
			}
		}
		else
		{
			for(int h = 0; h < imHP.hpMax; ++h)
			{
				final int ind = imHP.hp[h];
				final float val1 = (img1[ind]-avg1);
				final float val2 = (img2[ind]-avg2);
				std1 += val1*val1;
				std2 += val2*val2;
				mix += val1*val2;
			}
		}
		
		if( (std1 > 0) && (std2 > 0) )
		{
			return mix/(Math.sqrt(std1*std2));
		}
		else return 0.0; //neutral match to a constant image
	}
	
	
	void init_model(float[] model)
	{
		modelFloor = impAvg;
		Arrays.fill(model,0.0f);
		//check_floor(false);
	}
	
	//no bound checking here, do that before calling get_box
	float[] get_box(final float[] region, int rw, int rx, int ry, int bw, int bh)
	{
		final int barea = bw*bh;
		final int rh = region.length/rw;
		float[]  box= new float[barea];	
		
		HexPixels regHP = hex_pixels?(rw==impWidth ? imgHP : new HexPixels(rw/2)) : null;
		
		if(hex_pixels)
		{
			for(int bind = 0; bind < barea; ++bind)
			{
				final int bx = bind%bw;
				final int by = bind/bw;
				
				final int px = (rx+bx);
				final int py = (ry+by);
				
				//final int rind = px + py * rw;
				final int prind = regHP.periodic_qr_safe(px,py);
				box[bind] = region[prind];
			}	
		}
		else
		{
			for(int bind = 0; bind < barea; ++bind)
			{
				final int bx = bind%bw;
				final int by = bind/bw;
				
				final int px = (rx+bx+rw)%rw;
				final int py = (ry+by+rh)%rh;
				
				final int rind = px + py * rw;
				box[bind]=region[rind];
			}
		}
		return box;
	}
	
	//not hex_pixel aware and not used
	float[] get_diff_box(final float[] region2, final float[] region1, final int rw, final int rx, final int ry, final int bw, final int bh)
	{
		final int barea = bw*bh;
		final int rh = region1.length/rw;
		float[]  box= new float[barea];	
		for(int bind = 0; bind < barea; ++bind)
		{
			final int bx = bind%bw;
			final int by = bind/bw;
			
			final int px = (rx+bx+rw)%rw;
			final int py = (ry+by+rh)%rh;
			
			final int rind = px + py * rw;
			box[bind]=region2[rind]-region1[rind];
		}
		return box;
	}
	
	float[] get_diff(final float[] box1, final float[] box2)
	{
		final int area = box1.length;
		float[] diff = new float[area];	
		for(int ind = 0; ind < area; ++ind)
		{
			diff[ind]=box1[ind]-box2[ind];
		}
		return diff;
	}
	
	
	
	float fit_local_weight(float[] model, int pos, float weight)
	{
		//remember global state as we are going to work on a boxed copy
		//final double lastmodelFloor = modelFloor;
		final int lastnumAtoms = numAtoms;
		int bpos = pos;
		int bmw = impWidth;
		float[] box_model = null;
		float[] box_imp = imp_pix;
		float[] box_sim = sim_pix;
		if(!external_merit)
		{
			//crop a box from model and imp	
			int xm = pos % impWidth;
			int ym = pos / impWidth;
			
			if(hex_pixels)
			{
				//for positive as good as applying a bit mask with ^1
				xm = 2*(xm/2);
				ym = 2*(ym/2);
			}
			
			final int r = 1;
			
			int x1 = xm - r*(atomWidth) + (hex_pixels?1:0);
			int x2 = xm + r*(atomWidth);
			int y1 = ym - r*(atomHeight) + (hex_pixels?1:0);
			int y2 = ym + r*(atomHeight);
			
			
			if((!periodic) && (!hex_pixels))
			{
				if(x1 < 0) {x1 = 0;}
				if(y1 < 0) {y1 = 0;}
				if(x2 >= impWidth) {x2 = impWidth-1;}
				if(y2 >= impHeight) {y2 = impHeight-1;}
			}
			bmw = x2-x1+1;
			final int bmh = y2-y1+1;
			final int bmarea = bmw*bmh;
			
			final int bmx = xm-x1;
			final int bmy = ym-y1;
			bpos = bmx + bmy*bmw;
			//assume the atom with the weight is already there
			float best_weight = weight;
			box_model = get_box(model, impWidth, x1, y1, bmw, bmh);
			box_imp = get_box(imp_pix, impWidth, x1, y1, bmw, bmh);
			box_sim = sparse_floor_fold(box_model,bmw, modelFloor, a_pix, atomWidth);
		}
		else
		{
			box_model = Arrays.copyOf(model,model.length);
		}
		float best_weight = weight;
		double var0 = external_merit?model_merit:img_match(box_sim,box_imp,bmw);
		//float[] box_diff = get_diff(box_sim,box_imp);
		//double var0 = get_var(box_diff,0.0f);
		
		double f = 0.9; //hmm how far off may weights be?
		double fmax = 0.99;
		boolean shrinking = random.nextBoolean();
		int runs = -1;
		while( (f < fmax) && (++runs < 5) )
		{
			final float t_weight = (float)(shrinking ? f*best_weight : best_weight/f) ;
			enforce_atom(box_model, bmw, bpos, t_weight);
			box_sim = sparse_floor_fold(box_model,bmw, modelFloor, a_pix, atomWidth);
			//box_diff = get_diff(box_sim,box_imp);
			//double tmp_var = get_var(box_diff,0.0f);
			final double tmp_var = img_match(box_sim,box_imp,bmw);
			if(tmp_var > var0 )
			{
				best_weight = t_weight;
			}
			else
			{
				shrinking = !shrinking;
				f=Math.sqrt(f);
				//next or final put down will cancel the wrong weight anyways
			}	
		}
		
		if(external_merit)
		{	
			enforce_atom(box_model, bmw, bpos, 0.0f);
		}
		//modelFloor = lastmodelFloor;
		//numAtoms = lastnumAtoms;
		return (best_weight);
	
	}
	
	
	
	
	double estimate_improvement(float[] model, int pos, float weight)
	{
		//remember global state as we are going to work on a boxed copy
		//final double lastmodelFloor = modelFloor;
		//final int lastnumAtoms = numAtoms;
		int bmw = impWidth;
		int bpos = pos;
		float[] box_model = null;
		float[] box_sim = sim_pix;
		float[] box_imp = imp_pix;
		
		//crop a box from model and imp	
		if(!external_merit)
		{
			int xm = pos % impWidth;
			int ym = pos / impWidth;
			
			if(hex_pixels)
			{
				//for positive as good as applying a bit mask with ^1
				xm = 2*(xm/2);
				ym = 2*(ym/2);
			}
			
			final int r = 1;
			
			int x1 = xm - r*(atomWidth) + (hex_pixels?1:0);
			int x2 = xm + r*(atomWidth);
			int y1 = ym - r*(atomHeight) + (hex_pixels?1:0);
			int y2 = ym + r*(atomHeight);
			
			
			if( (!periodic) && (!hex_pixels) )
			{
				if(x1 < 0) {x1 = 0;}
				if(y1 < 0) {y1 = 0;}
				if(x2 >= impWidth) {x2 = impWidth-1;}
				if(y2 >= impHeight) {y2 = impHeight-1;}
			}
			bmw = x2-x1+1;
			final int bmh = y2-y1+1;
			final int bmarea = bmw*bmh;
			
			final int bmx = xm-x1;
			final int bmy = ym-y1;
			bpos = bmx + bmy*bmw;
			
			box_model = get_box(model, impWidth, x1, y1, bmw, bmh);
			box_sim = sparse_floor_fold(box_model,bmw, modelFloor, a_pix, atomWidth);
			box_imp = get_box(imp_pix, impWidth, x1, y1, bmw, bmh);
			
		}
		else
		{
			box_model = Arrays.copyOf(model,model.length);
		}
		//assume the area has been voided
		if(box_model[bpos] != 0.0)
		{
			throw new RuntimeException("Ups, estimate improvement got not voided area to work on");
		}
		
		final double void_merit = external_merit?model_merit:img_match(box_sim,box_imp,bmw); 
		//put down a the atom with the given weight //we could reuse box_sim here
		if(protect_atoms)
		{
			enforce_atom(box_model, bmw, bpos, weight);
		}
		else
		{
			put_down_atom(box_model, bmw, bpos, weight);
		}
		box_sim = sparse_floor_fold(box_model,bmw, modelFloor, a_pix, atomWidth);
		//double weight_match = img_match(weight_sim,box_imp,bmw);
		final double weight_merit = img_match(box_sim,box_imp,bmw); 
		
		//restore global state
		//modelFloor = lastmodelFloor;
		//numAtoms = lastnumAtoms;
		return (weight_merit - void_merit);
	}
	
	double estimate_weight(float[] model, int pos)
	{
		/*
		if(Float.isNaN(imp_pix[pos]))
		{	throw new RuntimeException(Macro.MACRO_CANCELED);}
		*/
		int bpos = pos;
		int bmw = impWidth;
		int bmh = impHeight;
		int bmx = pos%bmw;
		int bmy = pos/bmw;
		float[] box_model = null;
		float[] box_imp = imp_pix;
		float[] box_sim = sim_pix;
		
		if(!external_merit)
		{
			int xm = pos % impWidth;
			int ym = pos / impWidth;
			
			if(hex_pixels)
			{
				//for positive as good as applying a bit mask with ^1
				xm = 2*(xm/2);
				ym = 2*(ym/2);
			}
			
			final int r = 1;
			
			int x1 = xm - r*(atomWidth) + (hex_pixels?1:0);
			int x2 = xm + r*(atomWidth);
			int y1 = ym - r*(atomHeight) + (hex_pixels?1:0);
			int y2 = ym + r*(atomHeight);
			
			if((!periodic) && (!hex_pixels))
			{
				if(x1 < 0) {x1 = 0;}
				if(y1 < 0) {y1 = 0;}
				if(x2 >= impWidth) {x2 = impWidth-1;}
				if(y2 >= impHeight) {y2 = impHeight-1;}
			}
			
			bmw = x2-x1+1;
			bmh = y2-y1+1;
			final int bmarea = bmw*bmh;
			
			bmx = xm-x1;
			bmy = ym-y1;
			bpos = bmx + bmy*bmw;
			
			box_model = get_box(model, impWidth, x1, y1, bmw, bmh);
			box_imp = get_box(imp_pix, impWidth, x1, y1, bmw, bmh);
			box_sim = sparse_floor_fold(box_model,bmw, modelFloor, a_pix, atomWidth);
		}
		else
		{
			box_model = Arrays.copyOf(model,model.length);
		}
		
		double last_match = external_merit?model_merit:img_match(box_sim,box_imp,bmw); //the local match might differ from the global one
		
		double repelled = put_down_atom(box_model, bmw, bpos, 0.0 );
		if( ((!protect_atoms) && (deltaAtoms == -2)) || (protect_atoms && (deltaAtoms < -1) )) 
		{   //we probably tried to squeze in on top of a bond
			//modelFloor = lastmodelFloor;
			//numAtoms = lastnumAtoms;
			return Double.NaN; //dont suggest any action here
		}
		
		
		
		/* //This is not such a good idea on huge images with many atoms
		if((repelled > 0.0) && (repelled < avg_weight * peakWeightMin) ) 
		{
			//modelFloor = lastmodelFloor;
			//numAtoms = lastnumAtoms;
			return 0.0; //tell to eliminate an encountered weak atom
		}
		*/
		box_sim = sparse_floor_fold(box_model,bmw, modelFloor, a_pix, atomWidth);
		
		//determine the tentative weight here from simulation with voided atom
		//double guessed_weight = (protect_atoms ? repelled : 0.0);
		double guessed_weight = 0.0;
		
		double norm_inside = 0.0;
		//int pts_inside = 0;
		final int aarea = atomWidth * atomHeight; 
		final int maxpos = (atomHeight/2)*atomWidth + atomHeight/2;//should work for regular and hexpixels 
		final float aval_max = a_pix[maxpos];
		final float aval_min = 0.25f * aval_max;
		
		for(int aind = 0; aind < aarea; ++aind)
		{
			final float aval = a_pix[aind];
			final int ax = aind % atomWidth - atomWidth/2;
			final int ay = aind / atomWidth - atomHeight/2;
			
			final int xb = bmx + ax;
			final int yb = bmy + ay;
			
			if(	(xb >=0) && (yb >=0) &&
				(xb < bmw) && (yb < bmh) )
			{
				final int bind = xb+yb*bmw;
				final float dif = box_imp[bind]-box_sim[bind];
				if(aval >= aval_min)
				{
					guessed_weight += dif*aval;
					norm_inside += aval;
					//++pts_inside;
				}
			}	
		
		}
		guessed_weight *= (1.0/(aval_max*norm_inside)); //squaring a 2D Gauss should shrink sigma to half! 
		boolean reject_guess = false;
		reject_guess |= (repelled * guessed_weight < 0.0); //dont flip holes and atoms just because of inaccurate weights
		reject_guess |= ((!enable_holes) && (guessed_weight < 0.0) );
		reject_guess |= ( (guessed_weight > 0.0) && (guessed_weight < peakWeightMin * avg_peak_weight) && (numAtoms > 0) ); 
		reject_guess |= ( ( guessed_weight < 0.0) && (guessed_weight > peakWeightMin * avg_hole_weight) && (numHoles > 0) );
		
		if(reject_guess && protect_atoms)
		{	return Double.NaN;}
		
		double void_match = img_match(box_sim,box_imp,bmw);
		//if(!optimize_weights && (!((protect_atoms) && (repelled!=0.0))) )
		//{	guessed_weight = (guessed_weight > 0) ? avg_peak_weight : avg_hole_weight;}
		
		if(!reject_guess)
		{
			put_down_atom(box_model, bmw, bpos, guessed_weight );
			box_sim = sparse_floor_fold(box_model,bmw, modelFloor, a_pix, atomWidth);
		}
		double guessed_match = reject_guess?-1.0:img_match(box_sim,box_imp,bmw);
		
		if( (last_match > void_match) && (last_match > guessed_match) )//we tried to "move" atom to worse position
		{
			//modelFloor = lastmodelFloor;
			//numAtoms = lastnumAtoms;
			return Double.NaN; //no action
		}
		
		if( (void_match > last_match) && (void_match > guessed_match) && (!protect_atoms) ) // simply clear the area, maybe add an atom later
		{
			//modelFloor = lastmodelFloor;
			//numAtoms = lastnumAtoms;
			//return (repelled!=0.0) ? 0.0 : -1.0; //suggest to clear the area if there was something
			return 0.0;
		}
		
		if( (guessed_match > last_match) && (guessed_match > void_match) ) //suggest to nudge around that atom with full area matches
		{		
			//modelFloor = lastmodelFloor;
			//numAtoms = lastnumAtoms;
			return guessed_weight; 
		}
		//by default just try another spot
		//modelFloor = lastmodelFloor;
		//numAtoms = lastnumAtoms;
		return Double.NaN;
	}
	
	float enforce_atom(final float[] model,final int mW, final int pos, final double weight)
	{
		deltaAtoms=0;
		HexPixels imHP = null;
		if(hex_pixels)
		{
			if(mW == impWidth)
			{	imHP = imgHP;}
			else if(mW == 2*atomWidth)
			{	imHP = boxHP;}
			else //total overkill but should be cold code anyways
			{
				imHP = new HexPixels(mW/2);
			}
		}
		final boolean actual = (model == enh_pix);
		final float lifted = model[pos];
		model[pos] = (float)weight;
		
		if( (lifted > 0.0f) || (lifted < 0.0f) )
		{	--deltaAtoms;}
		if(weight > 0.0 || weight < 0.f)
		{	++deltaAtoms;}
		
		if( actual ) 
		{
			if(lifted > 0.0f)
			{	
				--numAtoms;
				total_peak_weight -= lifted;	
			}
			else if(lifted < 0.0f)
			{	
				--numHoles;
				total_hole_weight -= lifted;	
			}
			
			if(weight > 0.0)
			{	
				++numAtoms;
				total_peak_weight += weight;
			}
			else if(weight < 0.0)
			{	
				++numHoles;
				total_hole_weight += weight;
			}
			
			final int area = (hex_pixels ? imHP.hpMax : validArea);
			modelFloor = (impTotal-(total_peak_weight+total_hole_weight))/area;
			avg_peak_weight = (numAtoms > 0)?total_peak_weight/numAtoms:0.0;
			avg_hole_weight = (numHoles > 0)?total_hole_weight/numHoles:0.0;
			
			
		}
		//check_floor(false);
		return lifted;
	}
	
	
	
	
	//weight NaN triggers condensation
	float put_down_atom(final float[] model,final int mW, int pos, double weight) //weight = 0.0 for clearing
	{
		
		deltaAtoms = 0;
		HexPixels imHP = null;
		if(hex_pixels)
		{
			if(mW == impWidth)
			{	imHP = imgHP;}
			else if(mW == 2*atomWidth)
			{	imHP = boxHP;}
			else 
			{
				imHP = new HexPixels(mW/2);
			}
		}
		final boolean actual = (model == enh_pix);
		double lifted = 0.0;
		float biggest_present = 0.0f;
		int biggest_pos = -1; 
		
		if(hex_pixels)
		{
			
			
			for(int h = 0; h < tinyHP.hpMax; ++h)
			{
				final int qr = tinyHP.hp[h];
				final int[] txyz = tinyHP.to_xyz(qr);
				final int ind = imHP.shift_qr(pos,txyz);
				final float mval = model[ind];
				if( (mval > 0.0f) || (mval < 0.0f) )
				{
					if(actual)
					{	
						if( (mval > 0.0f)  )
						{	
							--numAtoms;
							total_peak_weight -= mval;
						}
						else if((mval < 0.0f) )
						{	
							--numHoles;
							total_hole_weight -= mval;
						}
						
					}
					 
					--deltaAtoms;
					lifted += mval;
					model[ind] = 0.0f;
					if(Math.abs(mval) > Math.abs(biggest_present))
					{
						biggest_present = mval;
						biggest_pos = ind;
					} 	
							
				}
			}
		
		}
		else // !hex_pixels
		{
			final int marea = model.length;
			final int mH = marea/mW;
			final int mx = pos%mW;
			final int my = pos/mW;
			
			final int rxr = (int) (atomSolidD * atomSolidD);
			final int ax1 = (int)(-atomSolidD-0.5);
			final int ax2 = (int)( atomSolidD+1);
			final int ay1 = (int)(-atomSolidD-0.5);
			final int ay2 = (int)( atomSolidD+1);
			
			
			for(int ay = ay1; ay < ay2; ++ay)
			{
				for(int ax = ax1; ax < ax2; ++ax)
				{	
					int px = mx + ax;
					int py = my + ay;
					
					if(periodic)
					{
						px = (px+mW) % mW;
						py = (py+mH) % mH;
					}
					
					if(	periodic || 
						( (px >= 0) && (px < mW) &&
						 (py >= 0) && (py < mH) )      )
					{
						final int ind = px + py * mW;
						final float mval = model[ind]; 
						if( ( (mval > 0.0f) || (mval < 0.0f) ) && (ax*ax + ay*ay <= rxr) )
						{
							if(actual)
							{	
								if( (mval > 0.0f) )
								{
									--numAtoms;
									total_peak_weight -= mval;
								}
								else if( (mval < 0.0f) )
								{
									--numHoles;
									total_hole_weight -= mval;	
								}
							}
							
							--deltaAtoms;
							lifted += mval;
							if(Math.abs(mval) > Math.abs(biggest_present))
							{
								biggest_present = mval;
								biggest_pos = ind;
							} 
							model[ind] = 0.0f;
							
						}
					}	
				}
			}
		}
		if(Double.isNaN(weight))
		{
			weight = lifted;
			pos = biggest_pos;
			//would be nice if we could only perform condensation at
			//the initial pos, but then we need an undo for the clearing out 
			
		}
		if( (pos != -1) && (weight != 0.0) )
		{
			model[pos] = (float)weight;
			if(weight != 0.0)
			{
				++deltaAtoms;
				if(actual)
				{	
					if(weight > 0.0)
					{	
						++numAtoms;
						total_peak_weight += weight;
						
					}
					else if (weight < 0.0)
					{	
						++numHoles;
						total_hole_weight += weight;
						
					}
				}
				
			}
		}
		if(actual)
		{
			final int area = (hex_pixels ? imgHP.hpMax : validArea);
			modelFloor = (impTotal-(total_peak_weight+total_hole_weight))/area;
			avg_peak_weight = (numAtoms > 0)?total_peak_weight/numAtoms:0.0;
			avg_hole_weight = (numHoles > 0)?total_hole_weight/numHoles:0.0;
		}
		
		//check_floor(false);
		return biggest_present;
	}
	
	void write_top()
	{
		try
		{
			FileInfo fi = imp.getOriginalFileInfo();
			String path = (fi!=null)?fi.directory:IJ.getDirectory("Choose a directory for "+ imp.getShortTitle() +".top file");
			File topfile = new File(path+imp.getShortTitle() + ".top");
			FileWriter fw = new FileWriter(topfile, false);
			PrintWriter pw = new PrintWriter(fw);
			
			pw.println("#creator: " + getClass().getSimpleName() + " source: " + impT);
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
				if( (bondi.left_ring == null) || (bondi.right_ring == null) ||
				    (!bondi.left_ring.is_interior) || (!bondi.right_ring.is_interior) )
				{	continue;}
				
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
				pw.println("#view of slice:" + (v+1) + " in " + impT);
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
			IJ.log("rewrote " + topfile.getName());
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	
	
	void write_sig(SurfaceMesh mmesh)
	{
		try
		{
			FileInfo fi = imp.getOriginalFileInfo();
			String path = (fi!=null)?fi.directory:IJ.getDirectory("Choose a directory for "+ imp.getShortTitle() +".sig file");
			File sigfile = new File(path+imp.getShortTitle() + ".sig");
			FileWriter fw = new FileWriter(sigfile, false);
			PrintWriter pw = new PrintWriter(fw);
			
			
			pw.println("#creator: " + getClass().getSimpleName() + " source: " + impT + " slice: " + impSlice);
			pw.println("#The signature lines are the polygons traversed when spiraling out from the choosen atom\n" +
					   "#the first line is always 5-6-7 each further line is another ring, the first ring in a \n" +
					    "#new line touches the last ring of the previous ring and shares at least an atom with the first ring");
			pw.println("ATOMS\t" + mmesh.fatoms.length);
			pw.println("BONDS\t" + mmesh.fbonds.length);
			pw.println("RINGS\t" + mmesh.frings.length);
			String[] sig_lines = mmesh.signature.sig_string.split("\\n");
			pw.println("SIGNATURE\t" + sig_lines.length);
			for(int i = 0; i < sig_lines.length; ++i)
			{
				pw.println(sig_lines[i]);  
			}
			pw.close();
			IJ.log("rewrote " + sigfile.getName());
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	
	String read_sig()
	{
		try
		{
			FileInfo fi = imp.getOriginalFileInfo();
			String path = (fi!=null)?fi.directory:IJ.getDirectory("Choose a directory for sig file");
			File sigfile = new File(path+imp.getShortTitle() + ".sig");
			
			IJ.log("reading " + sigfile.getPath());
			
			if(sigfile.isFile() && sigfile.setReadable(true) )
			{
				StringBuilder sig_builder = new StringBuilder(256);
				Scanner fs = new Scanner(new FileReader(sigfile));
				boolean sigcontent = false;
				while( fs.hasNextLine() )
				{
					String line = fs.nextLine();
					if(!line.equals("") && !line.startsWith("#"))
					{	
						
						boolean sigstart = (!sigcontent && line.startsWith("SIGNATURE"));
						if(sigstart)
						{
							sigcontent = true;
							continue;
						}
						
						if(sigcontent)
						{
							sig_builder.append(line);
							sig_builder.append('\n');
						}
					}
				}
				fs.close();
				return sig_builder.toString();
			}
			return null;
		}
		catch (Exception e)
		{
			System.out.println("Error: could not read signature from file");
			e.printStackTrace();
			return null;
		}	
	}
	
	void check_floor(boolean detailed)
	{
		double sum_imp = 0.0;
		double sum_enh = 0.0;
		double sum_sim = 0.0;
		int num_atoms = 0;
		int num_holes = 0;
		
		if(enh_pix == null || imp_pix == null)
		{	return;}
		
		if(hex_pixels)
		{
			for(int h = 0; h < imgHP.hpMax; ++h)
			{
				final int i = imgHP.hp[h];
				sum_enh += enh_pix[i];
				sum_imp += imp_pix[i];
				if(enh_pix[i] > 0.0f)
				{	++num_atoms;}
				if(enh_pix[i] < 0.0f)
				{	++num_holes;}
			}
			if(inverted)
			{
				sum_enh = -sum_enh;
				sum_imp = -sum_imp;
			}				
			sum_sim = sum_enh + modelFloor*imgHP.hpMax;
			boolean failed = false;
			if(num_atoms != numAtoms)
			{	
				System.out.println("num_atoms: " + num_atoms + " != numAtoms: " + numAtoms);	
				failed = true;	
			}
			if(num_holes != numHoles)
			
			{	
				System.out.println("num_holes: " + num_holes + " != numHoles: " + numHoles);	
				failed = true;
			}
			if(Math.abs((sum_sim-sum_imp)/(sum_sim+sum_imp)) > 0.0001 )
			{
				failed = true;
				
			}
			if(failed)
			{
				System.out.println("sum_sim: " + sum_sim + " != " + " sum_imp: " + sum_imp + "  ( sum_enh: " + sum_enh +
				 " modelFloor: " + modelFloor + " )");
			}
			else if(detailed)
			{	
				System.out.println("sum_sim: " + sum_sim + " ~= " + " sum_imp: " + sum_imp + "  ( sum_enh: " + sum_enh +
				 " modelFloor: " + modelFloor + " )");
			}	
			
			
			if(failed)
			{	throw new RuntimeException(Macro.MACRO_CANCELED);}
				
		}
		else
		{
			for(int i = 0; i < enh_pix.length; ++i)
			{
				sum_enh += enh_pix[i];
				if(!Float.isNaN(imp_pix[i]))
				{	sum_imp += imp_pix[i];}
				
				if(enh_pix[i] > 0.0f)
				{	++num_atoms;}
				else if(enh_pix[i] < 0.0f)
				{	++num_holes;}
			}	
			/*
			if(inverted)
			{
				sum_enh = sum_enh;	
				if(sum_imp > 0.0) //imp_pix are alreasy inverted by pull_pixels()
				{	sum_imp = -sum_imp;}	
			}
			*/	
			sum_sim = sum_enh + modelFloor*validArea;
			boolean failed = false;
			if(num_atoms != numAtoms)
			{	
				System.out.println("num_atoms: " + num_atoms + " != numAtoms: " + numAtoms);	
					failed = true;
			}
			
			if(num_holes != numHoles)
			{	
				System.out.println("num_holes: " + num_holes + " != numHoles: " + numHoles);	
				failed = true;
			}
			
			if(failed || Math.abs((sum_sim-sum_imp)/(sum_sim+sum_imp)) > 0.0001 )
			{
				System.out.println("sum_sim: " + sum_sim + " != " + " sum_imp: " + sum_imp + "  ( sum_enh: " + sum_enh + " modelFloor: " + modelFloor + " * validArea: " + validArea + " )");
				throw new RuntimeException(Macro.MACRO_CANCELED);
			}
			
		}
		return;
	}
		
}
