package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.LinkedHashSet;

import dna.Parser;
import dna.Timer;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextFile;
import fileIO.TextStreamWriter;

import align2.Shared;
import align2.Tools;

/**
 * Filters text lines by exact match or substring.
 * @author Brian Bushnell
 * @date Jul 6, 2015
 *
 */
public class FilterLines {

	public static void main(String[] args){
		Timer t=new Timer();
		FilterLines mb=new FilterLines(args);
		mb.process(t);
	}
	
	public FilterLines(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(a.equals("names")){
				if(b!=null){
					String[] x=b.split(",");
					for(String s : x){
						names.add(s);
					}
				}
			}else if(a.equals("lines") || a.equals("maxlines")){
				maxLines=Tools.parseKMG(b);
			}else if(a.equals("substrings") || a.equals("substring")){
				if(b==null){b="t";}
				if(b.equals("header")){
					lineSubstringOfName=true;
				}else if(b.equals("name")){
					nameSubstringOfLine=true;
				}else{
					nameSubstringOfLine=lineSubstringOfName=Tools.parseBoolean(b);
				}
			}else if(a.equals("prefix") || a.equals("prefixmode")){
				prefixMode=Tools.parseBoolean(b);
			}else if(a.equals("replace")){
				String[] split2=b.split(",");
				assert(split2.length==2);
				replace1=split2[0];
				replace2=split2[1];
			}else if(a.equals("casesensitive") || a.equals("case")){
				ignoreCase=!Tools.parseBoolean(b);
			}else if(a.equals("include") || a.equals("retain")){
				exclude=!Tools.parseBoolean(b);
			}else if(a.equals("exclude") || a.equals("remove")){
				exclude=Tools.parseBoolean(b);
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{
			String[] x=names.toArray(new String[names.size()]);
			names.clear();
			for(String s : x){
				Tools.addNames(s, names, true);
			}
		}
		if(ignoreCase){
			String[] x=names.toArray(new String[names.size()]);
			names.clear();
			for(String s : x){
				names.add(s.toLowerCase());
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println(out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, true);
	}
	
	void process(Timer t){
		
		
		final TextFile tf=new TextFile(ffin1);
		
		final TextStreamWriter tsw;
		if(out1!=null){
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
		}else{tsw=null;}
		
		long linesProcessed=0;
		
		long linesOut=0;
		long bytesOut=0;
		
		{
			for(String line0=tf.readLine(true); line0!=null; line0=tf.readLine(true)){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;

				String line=(ignoreCase ? line0.toLowerCase() : line0);
				if(replace1!=null){line=line.replace(replace1, replace2);}
				String prefix=null;
				if(prefixMode){
					for(int x=1; x<line.length(); x++){
						char c=line.charAt(x-1);
						char next=line.charAt(x);
						if(Character.isWhitespace(c)){
							prefix=line.substring(0, x).trim();
							break;
						}
					}
				}

				boolean keepThisLine;
				boolean match;
				{
					match=(names.contains(line) || (prefix!=null && names.contains(prefix)));
					if(!match && (nameSubstringOfLine || lineSubstringOfName)){
						for(String name : names){
							if((lineSubstringOfName && name.contains(line)) || (nameSubstringOfLine && line.contains(name))){match=true;}
							else if(prefix!=null && ((lineSubstringOfName && name.contains(prefix)) || (nameSubstringOfLine && prefix.contains(name)))){match=true;}
						}
					}
					keepThisLine=(match!=exclude);
				}
				
				//					assert(false) : names.contains(name)+", "+name+", "+prefix+", "+exclude;
				
				if(keepThisLine){
					tsw.println(line0);
					linesOut++;
					bytesOut+=line0.length();
				}
			}
		}

		errorState|=tf.close();
		errorState|=tsw.poisonAndWait();
		
		t.stop();
		
		double rpnano=linesProcessed/(double)(t.elapsed);
		
		outstream.println("\nTime:               "+t);
		outstream.println("Lines Processed:    "+linesProcessed+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Lines Out:          "+linesOut);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		assert(false) : "printOptions: TODO";
	}
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;

	private String out1=null;
	
	/*--------------------------------------------------------------*/
	
	private boolean exclude=true;
	private boolean nameSubstringOfLine=false;
	private boolean lineSubstringOfName=false;
	private boolean ignoreCase=true;
	private boolean prefixMode=false;
	private long maxLines=-1;

	private String replace1=null;
	private String replace2=null;
	
	private LinkedHashSet<String> names=new LinkedHashSet<String>();
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;

	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean useSharedHeader=false;
	
}
