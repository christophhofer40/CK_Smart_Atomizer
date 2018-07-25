package krami; 
import java.util.concurrent.atomic.AtomicInteger;

public class SharedMerit
{
	private double merit = Double.NaN;
	private AtomicInteger state = null;
	//-1 locked
	//0 ready and up to date
	//1 awaiting next merit
	//2 requesting update
	
	public SharedMerit()
	{
		state = new AtomicInteger(1); 	
	}

	public void promise_merit( )
	{
		try
		{
			while( state.get() == -1 )
			{	
				Thread.sleep(10);
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return;
		}
		
		state.set(1);
	}


	public double get_merit()
	{
		try
		{
			while( state.get() != 0 )
			{	Thread.sleep(10);}
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return Double.NaN;
		}
		return merit;
	}
	
	public void set_merit(double nmerit)
	{
		try
		{
			while( state.get() == -1 )
			{	Thread.sleep(10);}
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return;
		}
		
		if( state.get() == 1 )
		{	
			merit = nmerit;
			state.set(0);
		}
	}
	
	public int set_state(int nstate)
	{
		if(nstate == state.get())
		{	return nstate;}
		try
		{
			while( state.get() == -1 )
			{	
				if(nstate == -1)
				{	return -1;}
				Thread.sleep(10);
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return -2;
		}
		return state.getAndSet(nstate);
	}
	
	public int get_state()
	{
		return state.get();
	}
	


}
