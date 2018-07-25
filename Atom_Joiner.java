import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import java.awt.event.*;
import ij.plugin.filter.*;
import ij.plugin.filter.PlugInFilterRunner;

public class Atom_Joiner implements ExtendedPlugInFilter
{
	final int flags = (DOES_ALL | SNAPSHOT );
	static AtomJoinerCore atomjoiner = null;
	
	
	public int setup(final String arg, final ImagePlus imp) 
	{
		if(imp == null)
		{	
			IJ.noImage();
			return DONE;
		}
		else
		{
			atomjoiner = new AtomJoinerCore(flags);
			
			if(atomjoiner.init(imp)!=0)
			{	
				IJ.log("AtomJoiner: initialization failed");
				return DONE;
			}
			IJ.log("Draw Lines to connect/disconnect atoms");
			IJ.log("Draw Boxes to inspect connections");
			IJ.log("Draw Circles to edit atoms");
			return flags;
		}	 
	}
	
	public int showDialog(final ImagePlus imp, final String command, final PlugInFilterRunner pfr)
	{
		int flgs = atomjoiner.showDialog(imp, command, pfr);
		return flgs;
	}
	
	public void setNPasses (int nPasses) {}
	
	public void run(ImageProcessor ip) 
	{
		atomjoiner.run(ip);
		return;	
	}
	
	
	
}
