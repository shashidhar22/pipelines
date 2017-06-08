package jgi;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import kmer.KmerTableSet;

import stream.ByteBuilder;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadStreamWriter;

import dna.AminoAcid;
import dna.Data;
import dna.Parser;
import dna.Timer;
import fileIO.ReadWrite;
import fileIO.FileFormat;

import align2.ListNum;
import align2.LongList;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;
import align2.TrimRead;
import assemble.Tadpole;

/**
 * @author Brian Bushnell
 * @date Aug 14, 2012
 *
 */
public class BBMerge {
	
	
	public static void main(String[] args){
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
//		boolean old=Shared.USE_JNI;
//		Shared.USE_JNI=false; //TODO: This is for RQCFilter.  Can be removed.
		BBMerge mr=new BBMerge(args);
		mr.process();
//		Shared.USE_JNI=old;
		Read.VALIDATE_IN_CONSTRUCTOR=true;
	}
	
	
	private static void printOptions(){
		System.err.println("Please consult the shellscript for usage information.");
	}
	
	
	private static String[] preparse(String[] args){
		if(args==null){return new String[0];}
		int nulls=0;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);}
			
			
			if(a.equals("jni") || a.equals("usejni")){
				Shared.USE_JNI=Tools.parseBoolean(b);
			}else if(a.equals("showfullargs") || a.equalsIgnoreCase("showFullArgs")){
				showFullArgs=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("vstrict") || a.equals("verystrict")){
				vstrict=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("ustrict") || a.equals("ultrastrict")){
				ustrict=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("xstrict") || a.equals("hstrict") || a.equals("hyperstrict") || a.equals("maxstrict")){
				xstrict=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("strict")){
				strict=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("loose")){
				loose=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("vloose") || a.equals("veryloose")){
				vloose=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("uloose") || a.equals("ultraloose")){
				uloose=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("xloose") || a.equals("hloose") || a.equals("hyperloose") || a.equals("maxloose")){
				xloose=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("fast")){
				fast=Tools.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("default")){
				if(Tools.parseBoolean(b)){
					xstrict=ustrict=vstrict=strict=loose=vloose=uloose=xloose=fast=false;
				}
				args[i]=null;
				nulls++;
			}
		}
		
		if(nulls==0){return args;}
		ArrayList<String> args2=new ArrayList<String>(args.length-nulls+5);
		if(strict || vstrict || ustrict || xstrict){
			strict=true;
			loose=vloose=uloose=xloose=false;
			
			args2.add("maxbad=4");
			args2.add("margin=3");
			args2.add("minqo=8");
			args2.add("qualiters=2");
			
			if(xstrict){
				args2.add("ratiomode=t");
				args2.add("normalmode=t");
				args2.add("requireratiomatch=t");

				args2.add("minentropy=56");
				args2.add("minoverlap=14");
				args2.add("minoverlap0=3");
				
				args2.add("maxratio=0.055");
				args2.add("ratiomargin=12");
				args2.add("ratiooffset=0.65");
				args2.add("ratiominoverlapreduction=4");
				args2.add("efilter=2");
				args2.add("pfilter=0.25");
			}else if(ustrict){
				args2.add("ratiomode=t");
				args2.add("normalmode=t");
				args2.add("requireratiomatch=t");

				args2.add("minentropy=56");
				args2.add("minoverlap=14");
				args2.add("minoverlap0=3");
				
				args2.add("maxratio=0.045");
				args2.add("ratiomargin=12");
				args2.add("ratiooffset=0.5");
				args2.add("ratiominoverlapreduction=4");
				args2.add("efilter=2");
				args2.add("pfilter=0.03");
			}else if(vstrict){
				args2.add("ratiomode=t");
				args2.add("normalmode=f");

				args2.add("minentropy=52");
				args2.add("minoverlap=12");
				args2.add("minoverlap0=4");

				args2.add("maxratio=0.05");
				args2.add("ratiomargin=12");
				args2.add("ratiooffset=0.5");
				args2.add("ratiominoverlapreduction=4");
				args2.add("efilter=2");
				args2.add("pfilter=0.008");
			}else{
				args2.add("ratiomode=t");
				args2.add("normalmode=f");
				
				args2.add("minentropy=42");
				args2.add("minoverlap0=7");
				args2.add("minoverlap=11");
				
				args2.add("maxratio=0.075");
				args2.add("ratiomargin=7.5");
				args2.add("ratiooffset=0.55");
				args2.add("ratiominoverlapreduction=4");
				args2.add("efilter=4");
				args2.add("pfilter=0.0008");
			}
		}else if(loose || vloose || uloose || xloose){
			loose=true;
			strict=vstrict=ustrict=xstrict=false;
			args2.add("minoverlap=8");
			args2.add("minoverlap0=9");
			args2.add("qualiters=4");
			args2.add("mismatches=3");
			args2.add("margin=2");
			
			args2.add("ratiooffset=0.4");
			
			if(xloose){
				args2.add("owq=t");
				args2.add("ouq=t");
				args2.add("minentropy=22");
				args2.add("minoverlap=8");
				args2.add("minoverlap0=7");
				args2.add("maxratio=0.2");
				args2.add("mismatches=3");
				args2.add("ratiomargin=2");
				args2.add("normalmode=t");
				args2.add("pfilter=0.0000001");
				args2.add("efilter=8");
				args2.add("margin=2");
				args2.add("ratiominoverlapreduction=2");
			}else if(vloose || uloose){
				args2.add("owq=t");
				args2.add("ouq=t");
				if(uloose){
//					args2.add("maxratio=0.14");
//					args2.add("ratiomargin=2");
//					args2.add("normalmode=t");
//					args2.add("pfilter=0.0000001");
					
					
					args2.add("minoverlap=8");
					args2.add("minoverlap0=7");
					args2.add("mismatches=3");
					args2.add("margin=2");

					args2.add("ratiominoverlapreduction=2");
					args2.add("efilter=8");
					args2.add("maxratio=0.16");
					args2.add("ratiomargin=2.2");
					args2.add("pfilter=0.0000002");
					args2.add("minentropy=24");
				}else{
					args2.add("ratiominoverlapreduction=3");
					args2.add("maxratio=0.12");
					args2.add("ratiomargin=3");
					args2.add("pfilter=0.000004");
					args2.add("minentropy=28");
					args2.add("efilter=7.5");
					args2.add("ratiooffset=0.45");
				}
			}else{
				args2.add("maxratio=0.11");
				args2.add("ratiomargin=4.7");
				args2.add("ratiominoverlapreduction=2");
				args2.add("pfilter=0.00002");
				args2.add("efilter=8");
				args2.add("minentropy=30");
			}
		}else if(fast){
			args2.add("maxratio=0.08");
			args2.add("ratiomargin=2.5");
			args2.add("ratiominoverlapreduction=3");
			args2.add("pfilter=0.0002");
			args2.add("efilter=8");
			args2.add("minentropy=39");
			args2.add("mininsert0=50");
		}
		
		for(String s : args){
			if(s!=null){args2.add(s);}
		}
		return args2.toArray(new String[args2.size()]);
	}
	
	public BBMerge(String[] args){
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		System.err.println("BBMerge version "+version);
		
		{
			String[] args0=args;
			args=preparse(args);

			if(args0!=args && showFullArgs){
				System.err.println("Revised arguments: "+Arrays.toString(args)+"\n");
			}
		}
		
		Timer ttotal=new Timer();
		ttotal.start();
		
		in1=(args[0].indexOf('=')>0 ? null : args[0]);
		in2=(in1!=null && args.length>1 && args[1].indexOf('=')<0 ? args[1] : null);
		if(in2!=null && "null".equalsIgnoreCase(in2)){in2=null;}
		
		{
			if(in1!=null && !in1.contains(",") && !in1.startsWith("stdin.") && !in1.equals("stdin")){
				File f=new File(in1);
				if(!f.exists() || !f.isFile()){
					in1=null;
//					throw new RuntimeException(in1+" does not exist.");
				}
			}
			if(in2!=null && !in2.contains(",")){
				File f=new File(in2);
				if(!f.exists() || !f.isFile()){
					in2=null;
//					throw new RuntimeException(in2+" does not exist.");
				}else if(in1.equalsIgnoreCase(in2)){
					throw new RuntimeException("Both input files are the same.");
				}
			}
		}
		
		ReadWrite.MAX_ZIP_THREADS=Shared.threads()-1;
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		
		Shared.READ_BUFFER_LENGTH=Tools.max(Shared.READ_BUFFER_LENGTH, 400);
		
		boolean mm0set=false;
		
		Parser parser=new Parser();
		parser.trimq2=trimq;
		parser.minAvgQuality=minAvgQuality;
		parser.minReadLength=minReadLength;
		parser.maxReadLength=maxReadLength;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);}

			if(Parser.isJavaFlag(arg)){
				//jvm argument; do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseQualityAdjust(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("extra")){
				for(String s : b.split(",")){
					extra.add(s);
				}
			}else if(a.equals("useratio") || a.equals("ratio") || a.equals("ratiomode")){
				useRatioMode=Tools.parseBoolean(b);
			}else if(a.equals("usenormalmode") || a.equals("normalmode")){
				useNormalMode=Tools.parseBoolean(b);
			}else if(a.equals("requireratiomatch") || a.equals("rrm")){
				requireRatioMatch=Tools.parseBoolean(b);
			}else if(a.equals("maxratio")){
				MAX_RATIO=Float.parseFloat(b);
//				useRatioMode=true;
			}else if(a.equals("ratiomargin")){
				RATIO_MARGIN=Float.parseFloat(b);
//				useRatioMode=true;
			}else if(a.equals("ratiooffset")){
				RATIO_OFFSET=Float.parseFloat(b);
//				useRatioMode=true;
			}else if(a.equals("ratiominoverlapreduction")){
				MIN_OVERLAPPING_BASES_RATIO_REDUCTION=Integer.parseInt(b);
//				useRatioMode=true;
			}else if(a.equals("minentropy") || a.equals("entropy")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					minEntropyScore=Integer.parseInt(b);
				}else{
					useEntropy=Tools.parseBoolean(b);
				}
			}else if(a.equals("minoverlappingbases") || a.equals("minoverlapbases") || a.equals("minoverlap")){
				MIN_OVERLAPPING_BASES=Integer.parseInt(b);
			}else if(a.equals("minoverlappingbases0") || a.equals("minoverlapbases0") || a.equals("minoverlap0")){
				MIN_OVERLAPPING_BASES_0=Integer.parseInt(b);
			}else if(a.equals("minqo") || a.equals("minq")){
				MIN_QUALITY=(byte)Integer.parseInt(b);
			}else if(a.equals("maxq")){
				Read.MAX_MERGE_QUALITY=(byte)Integer.parseInt(b);
			}else if(a.equals("qualiters")){
				QUAL_ITERS=Tools.max(1, Integer.parseInt(b));
			}else if(a.equals("maxbadbases") || a.equals("maxbad") || a.equals("mismatches")){
				MAX_MISMATCHES=Integer.parseInt(b);
			}else if(a.equals("maxbadbases0") || a.equals("maxbad0") || a.equals("mismatches0")){
				MAX_MISMATCHES0=Integer.parseInt(b);
				mm0set=true;
			}else if(a.equals("margin")){
				MISMATCH_MARGIN=Integer.parseInt(b);
			}else if(a.equals("usemapping")){
				USE_MAPPING=Tools.parseBoolean(b);
			}else if(a.equals("bin")){
				bin=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=Shared.setThreads(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Tools.parseKMG(b);
			}else if(a.equals("outgood") || a.equals("outmerged") || a.equals("outm") || a.equals("out")){
				out1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outgood1") || a.equals("outmerged1") || a.equals("outm1") || a.equals("out1")){
				out1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outgood2") || a.equals("outmerged2") || a.equals("outm2") || a.equals("out2")){
				out2=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb") || a.equals("outu") || a.equals("outunmerged") || a.equals("outbad")){
				outb1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb1") || a.equals("outu1") || a.equals("outunmerged1") || a.equals("outbad1")){
				outb1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb2") || a.equals("outu2") || a.equals("outunmerged2") || a.equals("outbad2")){
				outb2=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outinsert") || a.startsWith("outi") || a.startsWith("outlength")){
				outinsert=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist") || a.equals("hist") || a.equals("histogram") || a.equals("ihist")){
				ihist=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outa") || a.equals("outadapter")){
				outAdapter=b;
				findAdapterSequence=(outAdapter!=null);
			}else if(a.equals("outc") || a.equals("outcardinality")){
				outCardinality=b;
//			}else if(a.equals("outputfailed")){
//				OUTPUT_FAILED=Tools.parseBoolean(b);outCardinality
			}else if(a.equals("mix")){
				MIX_BAD_AND_GOOD=Tools.parseBoolean(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				NONZERO_ONLY=Tools.parseBoolean(b);
			}else if(a.equals("showhiststats")){
				showHistStats=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				assert(false) : "verbose flag is static final; recompile to change it.";
//				verbose=Tools.parseBoolean(b);
			}else if(a.equals("join") || a.equals("merge")){
				join=Tools.parseBoolean(b);
				if(join){ecco=false;}
			}else if(a.equals("ecco") || a.equals("ecc") || a.equals("errorcorrect")){
				ecco=Tools.parseBoolean(b);
				if(ecco){join=false;}
			}else if(a.equals("tbo") || a.equals("trimbyoverlap")){
				trimByOverlap=Tools.parseBoolean(b);
			}else if(a.equals("useoverlap") || a.equals("usebases") || a.equals("matebyoverlap") || a.equals("matebybases")){
				MATE_BY_OVERLAP=Tools.parseBoolean(b);
			}
//			else if(a.startsWith("skipmated")){
//				SKIP_MATED_READS=Tools.parseBoolean(b);
//			}
			else if(a.equals("lowercase")){
				lowercaseAdapters=Tools.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("trimonfailure") || a.equals("tof")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					TRIM_ON_OVERLAP_FAILURE=Integer.parseInt(b);
				}else{
					TRIM_ON_OVERLAP_FAILURE=(Tools.parseBoolean(b) ? 1 : 0);
				}
			}else if(a.equals("overlapusingquality") || a.equals("ouq")){
				overlapUsingQuality=Tools.parseBoolean(b);
			}else if(a.equals("overlapwithoutquality") || a.equals("owoq") || a.equals("owuq") || a.equals("owq")){
				overlapWithoutQuality=Tools.parseBoolean(b);
			}else if(a.equals("maxExpectedErrors") || a.equals("mee") || a.equals("meefilter")){
				maxExpectedErrors=Float.parseFloat(b);
			}else if(a.equals("mi") || a.equals("minins") || a.equals("mininsert")){
				minInsert=Integer.parseInt(b);
			}else if(a.equals("mi0") || a.equals("mininsert0")){
				minInsert0=Integer.parseInt(b);
			}else if(a.equals("minprob")){
				minProb=Float.parseFloat(b);
				assert(minProb<1) : "minprob must be less than 1.  At 1, even kmers with 100% probablity of being error-free will be discarded.";
			}else if(a.equals("prealloc")){
				prealloc=Tools.parseBoolean(b);
			}else if(a.equals("prefilter")){
				prefilter=Tools.parseBoolean(b);
			}else if(a.equals("k")){
				kmerLength=Integer.parseInt(b);
			}else if(a.equals("efilter")){
				if(b==null || Character.isLetter(b.charAt(0))){
					boolean x=Tools.parseBoolean(b);
					if(!x){efilterRatio=0;}
				}else{
					efilterRatio=Float.parseFloat(b);
				}
				useEfilter=efilterRatio>0;
			}else if(a.equals("pfilter")){
				if(b==null || Character.isLetter(b.charAt(0))){
					boolean x=Tools.parseBoolean(b);
					if(!x){pfilterRatio=0;}
				}else{
					pfilterRatio=Float.parseFloat(b);
				}
			}else if(a.equals("efilteroffset")){
				efilterOffset=Float.parseFloat(b);
			}else if(a.equals("kfilter")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					filterCutoff=Integer.parseInt(b);
					useKFilter=filterCutoff>0;
				}else{
					useKFilter=Tools.parseBoolean(b);
				}
			}else if(a.equals("usequality")){
				useQuality=Tools.parseBoolean(b);
			}else if(a.equals("ignorequality")){
				useQuality=!Tools.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Tools.parseBoolean(b);
			}else if(a.equals("samplerate")){
				samplerate=Float.parseFloat(b);
				assert(samplerate<=1f && samplerate>=0f) : "samplerate="+samplerate+"; should be between 0 and 1";
			}else if(a.equals("sampleseed")){
				sampleseed=Long.parseLong(b);
			}else if(a.equals("recalibrate") || a.equals("recalibratequality") || a.equals("recal")){
				recalibrateQuality=Tools.parseBoolean(b);
			}else if(a.equals("recalpairnum") || a.equals("recalibratepairnum")){
				CalcTrueQuality.USE_PAIRNUM=Tools.parseBoolean(b);
			}else if(a.equals("path")){
				Data.setPath(b);
			}else if(a.equals("iupacton") || a.equals("itn")){
				iupacToN=Tools.parseBoolean(b);
			}
			
			//Extension parameters
			
			else if(a.equals("extendright") || a.equals("er") || a.equals("extend") || a.equals("extendright1") || a.equals("er1") || a.equals("extend1")){
				extendRight1=(int)Tools.parseKMG(b);
			}else if(a.equals("extendright2") || a.equals("er2") || a.equals("extend2")){
				extendRight2=(int)Tools.parseKMG(b);
			}else if(a.equals("extenditerations") || a.equals("iterations") || a.equals("ei") || a.equals("iters")){
				extendIterations=Tools.max(1, (int)Tools.parseKMG(b));
			}else if(a.equals("ecctadpole") || a.equals("ecct")){
				eccTadpole=Tools.parseBoolean(b);
			}else if(a.equals("shave") || a.equals("removedeadends")){
				shave=Tools.parseBoolean(b);
			}else if(a.equals("rinse") || a.equals("shampoo") || a.equals("removebubbles")){
				rinse=Tools.parseBoolean(b);
			}else if(a.equals("branchlower") || a.equals("branchlowerconst")){
				branchLowerConst=(int)Tools.parseKMG(b);
			}else if(a.equals("branchmult2")){
				branchMult2=(int)Tools.parseKMG(b);
			}else if(a.equals("branchmult1")){
				branchMult1=(int)Tools.parseKMG(b);
			}else if(a.equals("mincount") || a.equals("mincov") || a.equals("mindepth") || a.equals("min")){
				minCountSeed=minCountExtend=(int)Tools.parseKMG(b);
			}else if(a.equals("mindepthseed") || a.equals("mds") || a.equals("mincountseed") || a.equals("mcs")){
				minCountSeed=(int)Tools.parseKMG(b);
			}else if(a.equals("mindepthextend") || a.equals("mde") || a.equals("mincountextend") || a.equals("mce")){
				minCountExtend=(int)Tools.parseKMG(b);
			}else if(a.equals("ilb") || a.equals("ignoreleftbranches") || a.equals("ignoreleftjunctions") || a.equals("ibb") || a.equals("ignorebackbranches")){
				extendThroughLeftJunctions=Tools.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
//		assert(false) : ecco;
		minInsert=Tools.max(minInsert, MIN_OVERLAPPING_BASES);
		if(minInsert0<1){
			minInsert0=(Tools.max((int)(minInsert*0.75), 5, MIN_OVERLAPPING_BASES_0));
			int cap=(loose ? 50 : 35);
			minInsert0=Tools.min(cap, minInsert0);
		}
		minInsert0=Tools.min(minInsert, minInsert0);
		
		if(MATE_BY_OVERLAP && !useNormalMode && !useRatioMode){
			System.err.println("\n*** WARNING! Both normal and ratio mode were disabled; using normal mode. ***\n");
			useNormalMode=true;
		}
		
		loglog=(outCardinality==null ? null : new LogLog(1999, 8, 31, -1));
		
		{//Process parser fields
			Parser.processQuality();
			
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=(parser.trimq2!=null ? parser.trimq2 : new byte[] {parser.trimq});
			qtrim1=parser.qtrim1;
			qtrim2=(parser.qtrim2 || (parser.trimq2!=null && parser.trimq2.length>1));
			if(qtrim1==false && qtrim2==false){
				qtrim1=((qtrimLeft||qtrimRight)&&trimq[0]>=0);
			}
			minAvgQuality=parser.minAvgQuality;
			minAvgQualityBases=parser.minAvgQualityBases;
			minReadLength=Tools.max(1, parser.minReadLength);
			maxReadLength=(parser.maxReadLength<0 ? Integer.MAX_VALUE : parser.maxReadLength);
//			untrim=parser.untrim;
			
			forceTrimModulo=parser.forceTrimModulo;
			forceTrimLeft=parser.forceTrimLeft;
			forceTrimRight=parser.forceTrimRight;
			forceTrimRight2=parser.forceTrimRight2;
		}
		parseCustom=FASTQ.PARSE_CUSTOM;
		if(verbose){
//			assert(false) : "verbose flag is static final; recompile to change it.";
//			BBMergeOverlapper.verbose=true;
		}
		
		if(trimByOverlap){
			join=false;
		}
		
		if(!mm0set){
			MAX_MISMATCHES0=MAX_MISMATCHES+(loose ? 2 : 0);
		}
		
		if(MAX_MISMATCHES0<MAX_MISMATCHES){
			MAX_MISMATCHES0=MAX_MISMATCHES+(loose ? 2 : 0);
			System.err.println("MAX_MISMATCHES0 was set to "+MAX_MISMATCHES0+" to remain >=MAX_MISMATCHES");
		}
		
		if(MISMATCH_MARGIN>MAX_MISMATCHES){
			MISMATCH_MARGIN=MAX_MISMATCHES;
			System.err.println("MISMATCH_MARGIN was set to "+MISMATCH_MARGIN+" to remain >=MAX_MISMATCHES");
		}
		
		if(recalibrateQuality){CalcTrueQuality.initializeMatrices();}
		
		if(findAdapterSequence){
			for(int i=0; i<adapterCounts.length; i++){
				for(int j=0; j<adapterCounts[i].length; j++){
					adapterCounts[i][j]=new LongList(150);
				}
			}
		}
		
		if(in2==null && in1!=null && in1.contains("#") && !new File(in1).exists()){
			in2=in1.replaceFirst("#", "2");
			in1=in1.replaceFirst("#", "1");
		}
		
		if(out2==null && out1!=null && out1.contains("#")){
			out2=out1.replaceFirst("#", "2");
			out1=out1.replaceFirst("#", "1");
		}
		
		if(outb2==null && outb1!=null && outb1.contains("#")){
			outb2=outb1.replaceFirst("#", "2");
			outb1=outb1.replaceFirst("#", "1");
		}
		
		if(extendRight1>0 || extendRight2>0 || useKFilter || eccTadpole){
			ArrayList<String> list=new ArrayList<String>();
			list.add("in1="+in1);
			list.add("in2="+in2);
			if(extra.size()>0){
				StringBuilder sb=new StringBuilder("in=");
				String comma="";
				for(String s : extra){
					sb.append(comma);
					sb.append(s);
					comma=",";
				}
				list.add(sb.toString());
			}
			list.add("branchlower="+branchLowerConst);
			list.add("branchmult1="+branchMult1);
			list.add("branchmult2="+branchMult2);
			list.add("mincountseed="+minCountSeed);
			list.add("mincountextend="+minCountExtend);
			list.add("minprob="+minProb);
			list.add("k="+kmerLength);
			list.add("prealloc="+prealloc);
			list.add("prefilter="+prefilter);
			tadpole=Tadpole.makeTadpole(list.toArray(new String[0]), false);
		}else{
			tadpole=null;
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outb1, outb2, outinsert, ihist, outCardinality, outAdapter)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					out1+", "+out2+", "+outb1+", "+outb2+", "+outinsert+", "+ihist+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outb1, outb2, outinsert, ihist)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		if(in2!=null){
			assert(!in1.equalsIgnoreCase(in2));
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
		}else{
			FASTQ.TEST_INTERLEAVED=true;
			FASTQ.FORCE_INTERLEAVED=true;
		}
		
		if(THREADS<1){THREADS=Shared.threads();}
		
		useMEEfilter=maxExpectedErrors>0;
		
		Read.VALIDATE_IN_CONSTRUCTOR=(THREADS<16);
	}
	
	void process(){
		Timer ttotal=new Timer();
		ttotal.start();
		
		if(tadpole!=null){
			Timer tload=new Timer();
			Tadpole.showSpeed=false;
			KmerTableSet.showSpeed=false;
			long kmers=tadpole.loadKmers(tload);
			tload.stop();
			System.err.println();
			
			if(shave || rinse){
				tload.start();
				long removed=tadpole.shaveAndRinse(tload, shave, rinse, true);
				tload.stop();
				System.err.println();
			}
			
//			System.err.println("Loaded "+kmers+" kmers in "+tload);
		}
		
		runPhase(join, maxReads, false);
		
		double stdev=0;
		if(histTotal!=null){
			stdev=Tools.standardDeviationHistogram(histTotal);
		}
		
		final long sum=correctCountTotal+incorrectCountTotal;
		final double divp=100d/readsProcessedTotal;
		final double div2=100d/sum;
		
		writeHistogram(ihist, sum*divp);
		
		if(outAdapter!=null){
			assert(findAdapterSequence);
			writeAdapterConsensus(outAdapter, adapterCounts);
		}
		
		if(outCardinality!=null){
			ReadWrite.writeString(loglog.cardinality()+"\n", outCardinality);
		}
		
		ttotal.stop();
		System.err.println("Total time: "+ttotal+"\n");
		
		System.err.println("Pairs:               \t"+readsProcessedTotal);
		System.err.println("Joined:              \t"+sum+String.format((sum<10000 ? "       " : "   ")+"\t%.3f%%", sum*divp));
		System.err.println("Ambiguous:           \t"+ambiguousCountTotal+String.format((ambiguousCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", ambiguousCountTotal*divp));
		System.err.println("No Solution:         \t"+noSolutionCountTotal+String.format((noSolutionCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", noSolutionCountTotal*divp));
		if(minInsert>0){System.err.println("Too Short:           \t"+tooShortCountTotal+String.format((tooShortCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", tooShortCountTotal*divp));}
		if(maxReadLength<Integer.MAX_VALUE){System.err.println("Too Long:            \t"+tooLongCountTotal+String.format((tooLongCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", tooLongCountTotal*divp));}
		
		if(extendRight1>0 || extendRight2>0){
			double dive=100d/extensionsAttempted;
			System.err.println("Fully Extended:      \t"+fullyExtendedTotal+String.format((fullyExtendedTotal<10000 ? "       " : "   ")+"\t%.3f%%", fullyExtendedTotal*dive));
			System.err.println("Partly Extended:     \t"+partlyExtendedTotal+String.format((partlyExtendedTotal<10000 ? "       " : "   ")+"\t%.3f%%", partlyExtendedTotal*dive));
			System.err.println("Not Extended:        \t"+notExtendedTotal+String.format((notExtendedTotal<10000 ? "       " : "   ")+"\t%.3f%%", notExtendedTotal*dive));
		}
		
		if(parseCustom){
			System.err.println();
			System.err.println("Correct:             \t"+correctCountTotal+String.format((correctCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", correctCountTotal*divp)+String.format("   \t%.3f%% of merged", correctCountTotal*div2));
			System.err.println("Incorrect:           \t"+incorrectCountTotal+String.format((incorrectCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", incorrectCountTotal*divp)+String.format("   \t%.3f%% of merged", incorrectCountTotal*div2));
			double snr=Tools.max(correctCountTotal, 0.001)/(Tools.max(incorrectCountTotal, 0.001));
			double snrDB=Tools.mid(-20, 80, 10*Math.log10(snr));
			System.err.println("SNR:                 \t"+String.format("%.3f dB", snrDB));
			System.err.println();
			System.err.println("Avg Insert Correct:  \t"+String.format("%.1f", (insertSumCorrectTotal)*1d/(correctCountTotal)));
			System.err.println("Avg Insert Incorrect:\t"+String.format("%.1f", (insertSumIncorrectTotal)*1d/(incorrectCountTotal)));
		}
		
		System.err.println("\nAvg Insert:          \t"+String.format("%.1f", (insertSumCorrectTotal+insertSumIncorrectTotal)*1d/(correctCountTotal+incorrectCountTotal)));
		System.err.println("Standard Deviation:  \t"+String.format("%.1f", stdev));
		System.err.println("Mode:                \t"+Tools.calcMode(histTotal));
		
		System.err.println();
		System.err.println("Insert range:        \t"+insertMinTotal+" - "+insertMaxTotal);
		System.err.println("90th percentile:     \t"+Tools.percentile(histTotal, .9));
		System.err.println("75th percentile:     \t"+Tools.percentile(histTotal, .75));
		System.err.println("50th percentile:     \t"+Tools.percentile(histTotal, .5));
		System.err.println("25th percentile:     \t"+Tools.percentile(histTotal, .25));
		System.err.println("10th percentile:     \t"+Tools.percentile(histTotal, .1));
	}
	
	public static void writeHistogram(String fname, double percentMerged){
		if(fname==null){return;}
		StringBuilder sb=new StringBuilder();

		if(showHistStats){
			sb.append("#Mean\t"+String.format("%.3f", Tools.averageHistogram(histTotal))+"\n");
			sb.append("#Median\t"+Tools.percentile(histTotal, 0.5)+"\n");
			sb.append("#Mode\t"+Tools.calcMode(histTotal)+"\n");
			sb.append("#STDev\t"+String.format("%.3f", Tools.standardDeviationHistogram(histTotal))+"\n");
			sb.append("#PercentOfPairs\t"+String.format("%.3f", percentMerged)+"\n");
		}
		sb.append("#InsertSize\tCount\n");
		for(int i=0; i<histTotal.length && i<=insertMaxTotal; i+=bin){
			int x=0;
			int y=0;
			for(int j=i; j<i+bin && j<histTotal.length; j++){
				x+=histTotal[j];
				y++;
			}
			x=(x+bin-1)/y;
			if(x>0 || !NONZERO_ONLY){
				sb.append(i+"\t"+x+"\n");
			}
		}
		ReadWrite.writeStringInThread(sb, fname);
	}
	
	public static void writeAdapterConsensus(String fname, LongList[][] matrix){
		StringBuilder sb=new StringBuilder();
		{
			sb.append(">Read1_adapter\n");
			StringBuilder adapter=new StringBuilder();
			LongList[] lists=matrix[0];
			long max=0;
			int lastBase=-1;
			for(int i=0; true; i++){
				long a=lists[0].get(i);
				long c=lists[1].get(i);
				long g=lists[2].get(i);
				long t=lists[3].get(i);
				long sum=(a+c+g+t);
				max=Tools.max(max, sum);
				if(sum==0 || (sum<10 && sum<=max/1000) || (max>100 && sum<8)){break;}
				long thresh=(max>100 ? 4+(sum*2)/3 : (sum*2)/3);
				if(a>thresh){
					adapter.append('A');
					lastBase=i;
				}else if(c>thresh){
					adapter.append('C');
					lastBase=i;
				}else if(g>thresh){
					adapter.append('G');
					lastBase=i;
				}else if(t>thresh){
					adapter.append('T');
					lastBase=i;
				}else{
					adapter.append('N');
				}
			}
			if(lastBase<0){sb.append('N');}
			else{
				for(int i=0; i<=lastBase; i++){
					sb.append(adapter.charAt(i));
				}
			}
			sb.append('\n');
		}
		if(matrix.length>1){
			sb.append(">Read2_adapter\n");
			StringBuilder adapter=new StringBuilder();
			LongList[] lists=matrix[1];
			long max=0;
			int lastBase=-1;
			for(int i=0; true; i++){
				long a=lists[0].get(i);
				long c=lists[1].get(i);
				long g=lists[2].get(i);
				long t=lists[3].get(i);
				long sum=(a+c+g+t);
				max=Tools.max(max, sum);
				if(sum==0 || (sum<10 && sum<=max/1000) || (max>100 && sum<8)){break;}
				long thresh=(max>100 ? 5+(sum*2)/3 : (sum*2)/3);
				if(a>thresh){
					adapter.append('A');
					lastBase=i;
				}else if(c>thresh){
					adapter.append('C');
					lastBase=i;
				}else if(g>thresh){
					adapter.append('G');
					lastBase=i;
				}else if(t>thresh){
					adapter.append('T');
					lastBase=i;
				}else{
					adapter.append('N');
				}
			}
			if(lastBase<0){sb.append('N');}
			else{
				for(int i=0; i<=lastBase; i++){
					sb.append(adapter.charAt(i));
				}
			}
			sb.append('\n');
		}
		ReadWrite.writeString(sb, fname);
	}
	
	public void runPhase(boolean join, long maxReads, boolean perfectonly){
		
		Timer talign=new Timer();
		
		ConcurrentReadOutputStream rosgood=null;
		ConcurrentReadOutputStream rosbad=null;
		ConcurrentReadOutputStream rosinsert=null;
		
		if(out1!=null){
			if(join==true){
				if(out2==null){System.err.println("Writing mergable reads merged.");}
				else{
					System.err.println("WARNING: 2 output files specified even though 'merge=true'.  out2 will be ignored.");
					out2=null;
				}
			}else{
				if(out2==null){System.err.println("Writing mergable reads interleaved.");}
				else{System.err.println("Writing mergable reads unmerged in two files.");}
			}
			
			final FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			final FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(16, 2*THREADS);
			rosgood=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			rosgood.start();
		}
		
		if(outb1!=null){

			final FileFormat ff1=FileFormat.testOutput(outb1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			final FileFormat ff2=FileFormat.testOutput(outb2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(16, 2*THREADS);
			rosbad=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			rosbad.start();
		}
		
		if(outinsert!=null){
			final int buff=Tools.max(16, 2*THREADS);
			
			String out1=outinsert.replaceFirst("#", "1");

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			
			ReadStreamWriter.HEADER=header();
			final FileFormat ff=FileFormat.testOutput(out1, FileFormat.ATTACHMENT, ".info", true, overwrite, append, ordered);
			rosinsert=ConcurrentReadOutputStream.getStream(ff, null, null, null, buff, null, false);
			rosinsert.start();
		}
		
		
		if(rosgood!=null || rosbad!=null || rosinsert!=null){
			System.err.println("Started output threads.");
		}
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.setSampleRate(samplerate, sampleseed);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
//		assert(paired);//Fails on empty files.
		if(verbose){System.err.println("Paired: "+paired);}
		
		talign.start();
		
		
		MateThread[] pta=new MateThread[THREADS];
		for(int i=0; i<pta.length; i++){
			pta[i]=new MateThread(cris, rosgood, rosbad, rosinsert, join, trimByOverlap);
			pta[i].start();
		}

		insertMinTotal=999999999;
		insertMaxTotal=0;
		
		readsProcessedTotal=0;
		matedCountTotal=0;
		correctCountTotal=0;
		ambiguousCountTotal=0;
		tooShortCountTotal=0;
		tooLongCountTotal=0;
		incorrectCountTotal=0;
		noSolutionCountTotal=0;
		insertSumCorrectTotal=0;
		insertSumIncorrectTotal=0;
		
		Arrays.fill(histTotal, 0);
		
		for(int i=0; i<pta.length; i++){
			MateThread ct=pta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
				readsProcessedTotal+=ct.pairsProcessed;
				matedCountTotal+=ct.matedCount;
				correctCountTotal+=ct.correctCount;
				ambiguousCountTotal+=ct.ambiguousCount;
				tooShortCountTotal+=ct.tooShortCount;
				tooLongCountTotal+=ct.tooLongCount;
				incorrectCountTotal+=ct.incorrectCount;
				noSolutionCountTotal+=ct.noSolutionCount;
				insertSumCorrectTotal+=ct.insertSumCorrect;
				insertSumIncorrectTotal+=ct.insertSumIncorrect;
				
				fullyExtendedTotal+=ct.fullyExtendedT;
				partlyExtendedTotal+=ct.partlyExtendedT;
				notExtendedTotal+=ct.notExtendedT;
				extensionsAttempted+=ct.extensionsAttemptedT;

				insertMinTotal=Tools.min(ct.insertMin, insertMinTotal);
				insertMaxTotal=Tools.max(ct.insertMax, insertMaxTotal);
				
//				System.err.println(ct.insertMin+", "+ct.insertMax);
				
				if(ct.hist!=null){
					for(int h=0; h<ct.hist.length; h++){
						histTotal[h]+=ct.hist[h];
					}
				}
				
				if(findAdapterSequence){
					LongList[][] adapterCountsT=ct.adapterCountsT;
					for(int x=0; x<adapterCounts.length; x++){
						for(int y=0; y<adapterCounts[x].length; y++){
							adapterCounts[x][y].add(adapterCountsT[x][y]);
						}
					}
				}
			}
		}
		
//		System.err.println("Finished reading");
		errorState|=ReadWrite.closeStreams(cris, rosgood, rosbad, rosinsert);
		
		talign.stop();
//		System.err.println("Align time: "+talign);
	}
	
	public static final float mergeableFraction(String fname1, String fname2, long numReads, float samplerate){
		long[] hist=makeInsertHistogram(fname1, fname2, numReads, samplerate);
		if(hist==null || hist.length<2){return 0;}
		long sum=Tools.sum(hist);
		return sum<1 ? 0 : (sum-hist[0])/(float)sum;
	}
	
	public static final long[] makeInsertHistogram(String fname1, String fname2, long numReads, float samplerate){
		assert(fname1!=null);
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(fname1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(fname2, FileFormat.FASTQ, null, true, true);
			if(ff1.stdio()){return null;}
			assert(!ff1.stdio()) : "Standard in is not allowed as input when calculating insert size distributions for files.";
			cris=ConcurrentReadInputStream.getReadInputStream(numReads, true, ff1, ff2);
			cris.setSampleRate(samplerate, 1);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
			if(!cris.paired()){
				ReadWrite.closeStreams(cris);
				return null;
			}
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert(r.mate!=null);
		}

		LongList ll=new LongList(500);
		while(reads!=null && reads.size()>0){

			for(Read r1 : reads){
				int x=findOverlapLoose(r1, r1.mate, false);
				if(x>0){ll.increment(x, 1);}
				else{ll.increment(0, 1);}
			}
			cris.returnList(ln.id, ln.list.isEmpty());
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln.id, ln.list.isEmpty());
		ReadWrite.closeStreams(cris);
		return ll.toArray();
	}

	/** Returns the insert size as calculated by overlap, or -1 */
	public static final int findOverlapStrict(final Read r1, final Read r2, boolean ecc){
		final float maxRatio=0.06f;
		final float ratioMargin=10f;
		final float ratioOffset=0.5f;
		
		final float efilterRatio=2f;
		final float efilterOffset=0.45f;
		final float pfilterRatio=0.008f;

		final int minOverlap=8;
		final int minOverlap0=4;
		final int minInsert=50;
		final int minInsert0=35;
		final int entropy=42;
		
		final int x=findOverlap(r1, r2, ecc,
				minOverlap, minOverlap0, minInsert, minInsert0, entropy,
				maxRatio, ratioMargin, ratioOffset,
				efilterRatio, efilterOffset, pfilterRatio);
		return x;
	}

	/** Returns the insert size as calculated by overlap, or -1 */
	public static final int findOverlapLoose(final Read r1, final Read r2, boolean ecc){
		
		final float maxRatio=0.12f;
		final float ratioMargin=3f;
		final float ratioOffset=0.45f;
		
		final float efilterRatio=7.5f;
		final float efilterOffset=0.55f;
		final float pfilterRatio=0.000004f;

		final int minOverlap=5;
		final int minOverlap0=6;
		final int minInsert=16;
		final int minInsert0=16;
		final int entropy=28;
		
		final int x=findOverlap(r1, r2, ecc,
				minOverlap, minOverlap0, minInsert, minInsert0, entropy,
				maxRatio, ratioMargin, ratioOffset,
				efilterRatio, efilterOffset, pfilterRatio);
		return x;
	}
	
	/** Returns the insert size as calculated by overlap, or -1 */
	public static final int findOverlap(final Read r1, final Read r2, final boolean ecc,
			int minOverlap, final int minOverlap0, final int minInsert, final int minInsert0, final int entropy,
			final float maxRatio, final float ratioMargin, final float ratioOffset,
			final float efilterRatio, final float efilterOffset, final float pfilterRatio){
		
		assert(r1!=null && r2!=null);
		if(!r1.validated()){r1.validate(true);}
		if(!r2.validated()){r2.validate(true);}

		final int len1=r1.length(), len2=r2.length();
		final int minlen=Tools.min(len1, len2);
		
		if(minlen<MIN_OVERLAPPING_BASES || minlen<minInsert){
			return -1;
		}
		
		int[] rvector=localRvector.get();
		if(rvector==null){
			rvector=new int[5];
			localRvector.set(rvector);
		}
		
		r2.reverseComplement();
		
		int bestInsert=-1;
		int bestBad=999999;
		boolean ambig, tooShort=false;
		
		if(USE_MAPPING && r1.chrom==r2.chrom && r1.start<r1.stop && r1.mapped() && r2.mapped()){
			bestBad=0;
			bestInsert=Read.insertSizeMapped(r1, r2, ignoreMappingStrand);
			ambig=false;
		}else{
			if(entropy>0){
				int a=BBMergeOverlapper.calcMinOverlapByEntropy(r1.bases, 3, null, entropy);
				int b=BBMergeOverlapper.calcMinOverlapByEntropy(r2.bases, 3, null, entropy);
				minOverlap=Tools.max(MIN_OVERLAPPING_BASES, Tools.max(a, b));
			}else{minOverlap=MIN_OVERLAPPING_BASES;}
			if(verbose){System.err.println("minOverlap: "+minOverlap);}
			
			rvector[4]=0;

			int x=BBMergeOverlapper.mateByOverlapRatio(r1, r2, null, null, rvector, minOverlap0, minOverlap, 
					minInsert0, minInsert, maxRatio, ratioMargin, ratioOffset, 0.95f, 0.95f, false);
			bestInsert=x;
			bestBad=rvector[2];
			ambig=(x>-1 ? rvector[4]==1 : false);
		}
		
		//TODO:  Crucial!  This line can vastly reduce merge rate, particularly if quality values are inaccurate. 
		if(bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null){
			float bestExpected=BBMergeOverlapper.expectedMismatches(r1, r2, bestInsert);
			if((bestExpected+efilterOffset)*efilterRatio<bestBad){ambig=true;}
			if(verbose){System.err.println("Result after efilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
		}
		
		//TODO:  Crucial!  This line can vastly reduce merge rate, particularly if quality values are inaccurate. 
		if(pfilterRatio>0 && bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null){
			float probability=BBMergeOverlapper.probability(r1, r2, bestInsert);
			if(probability<pfilterRatio){bestInsert=-1;}
			if(verbose){System.err.println("Result after pfilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
		}
		
		tooShort=(!ambig && bestInsert>0 && bestInsert<minInsert);
		
		if(ecc && bestInsert>-1 && !ambig && !tooShort){
			errorCorrectWithInsert(r1, r2, bestInsert);
		}
		
		if(r2!=null){r2.reverseComplement();}
		if(!ambig && bestInsert>-1){r1.setInsert(bestInsert);}
		
		return ambig ? -1 : bestInsert;
	}
	
	public static int errorCorrectWithInsert(Read r1, Read r2, int insert){
		assert(insert>0);
		int errors=0;
		Read joined=r1.joinRead(insert);
		
		if(joined!=null && joined.length()>0){
			final int lenj=joined.length();
			final int lim1=Tools.min(joined.length(), r1.length());
			final int lim2=lenj-Tools.min(r2.length(), lenj);

			r1.bases=Arrays.copyOfRange(joined.bases, 0, lim1);
			r1.quality=(r1.quality==null ? null : Arrays.copyOfRange(joined.quality, 0, lim1));

			r2.bases=Arrays.copyOfRange(joined.bases, lim2, lenj);
			r2.quality=(r2.quality==null ? null : Arrays.copyOfRange(joined.quality, lim2, lenj));
		}
		return errors;
	}

	public static String header(){
		return "#id\tnumericID\tinsert\tstatus\tmismatches\n";
	}
	
	private void qtrim(Read r1, Read r2, int iter){
		if(false /*untrim*/){
			TrimRead.trim(r1, qtrimLeft, qtrimRight, trimq[iter], 1);
			TrimRead.trim(r2, qtrimLeft, qtrimRight, trimq[iter], 1);
		}else{
			TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq[iter], 1);
			TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq[iter], 1);
		}
	}
	
	private class MateThread extends Thread{
		
		
		public MateThread(ConcurrentReadInputStream cris_, ConcurrentReadOutputStream rosgood_, ConcurrentReadOutputStream rosbad_, ConcurrentReadOutputStream rosi_,
				boolean joinReads_, boolean trimByOverlap_) {
			cris=cris_;
			rosgood=rosgood_;
			rosbad=rosbad_;
			rosi=rosi_;
			joinReads=joinReads_;
			trimReadsByOverlap=trimByOverlap_;
			
			if(useEntropy){
				kmerCounts=new short[1<<(2*entropyK)];
			}else{
				kmerCounts=null;
			}
			
			if(findAdapterSequence){
				for(int i=0; i<adapterCountsT.length; i++){
					for(int j=0; j<adapterCountsT[i].length; j++){
						adapterCountsT[i][j]=new LongList(150);
					}
				}
			}
		}
		
		
		@Override
		public void run(){
			processReads();
		}
		
		private void processReads() {
			assert(USE_MAPPING || MATE_BY_OVERLAP);
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(r.mate!=null);
			}
			
			final byte[][] originals=(rosbad!=null || (rosgood!=null && (!join || MIX_BAD_AND_GOOD))) ? new byte[4][] : null;
			while(reads!=null && reads.size()>0){
				
				ArrayList<Read> listg=(rosgood==null /*&& rosi==null*/ ? null : new ArrayList<Read>(reads.size()));
				ArrayList<Read> listb=(rosbad==null ? null : new ArrayList<Read>(reads.size()));
				
				if(loglog!=null){
					for(Read r1 : reads){loglog.hash(r1);}
				}
				
				for(Read r1 : reads){
					int bestInsert=findOverlapInThread(r1, originals, listg, listb);
				}
				
				if(rosgood!=null){rosgood.add(listg, ln.id);}
				if(rosi!=null){
					//This prints both merged and unmerged reads
					for(Read r1 : reads){//Legacy outinsert support
						StringBuilder sb=new StringBuilder(40);
						sb.append(r1.id==null ? r1.numericID+"" : r1.id).append('\t');
						sb.append(r1.numericID).append('\t');
						final int bestInsert=r1.insert();
						sb.append(bestInsert>=0 ? bestInsert : -1);
						sb.append('\t');

						if(bestInsert==RET_NO_SOLUTION){sb.append('F');}//Failed
						else if(bestInsert==RET_AMBIG){sb.append('A');} //Ambiguous
						else if(bestInsert==RET_SHORT){sb.append('S');} //Short
						else{
							if(r1.errors>0){sb.append('I');}//Imperfect
							else{sb.append('P');}//Perfect
							sb.append('\t');
							sb.append(r1.errors);
						}
						r1.obj=sb;
					}
					rosi.add(reads, ln.id);
				}
				if(rosbad!=null){rosbad.add(listb, ln.id);}
				
				//			System.err.println("returning list");
				cris.returnList(ln.id, ln.list.isEmpty());
				//			System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
				//			System.err.println("reads: "+(reads==null ? "null" : reads.size()));
			}
			cris.returnList(ln.id, ln.list.isEmpty());
		}
		
		private int findOverlapInThread(final Read r1, final byte[][] originals, ArrayList<Read> listg, ArrayList<Read> listb){

			final Read r2=r1.mate;
			final int trueSize=r1.insert();
			
			final int bestInsert=processReadPair(r1, r2);
			
			if(originals!=null){
				if(eccTadpole){
					originals[0]=r1.bases;
					originals[1]=r1.quality;
					originals[2]=r2.bases;
					originals[3]=r2.quality;
				}else{
					originals[0]=(r1.bases==null ? null : r1.bases.clone());
					originals[1]=(r1.quality==null ? null : r1.quality.clone());
					originals[2]=(r2.bases==null ? null : r2.bases.clone());
					originals[3]=(r2.quality==null ? null : r2.quality.clone());
				}
			}
			
			Read joined=null;
			
			if(bestInsert>0){
				if(bestInsert==trueSize){correctCount++;insertSumCorrect+=bestInsert;}
				else{incorrectCount++;insertSumIncorrect+=bestInsert;}
				r1.setInsert(bestInsert);
				insertMin=Tools.min(bestInsert, insertMin);
				insertMax=Tools.max(bestInsert, insertMax);
				hist[Tools.min(bestInsert, hist.length-1)]++;
				if(joinReads){
					r2.reverseComplement();
					joined=r1.joinRead(bestInsert);
					r2.reverseComplement();
					assert(joined.length()==bestInsert);
				}else if(ecco){
					r2.reverseComplement();
					errorCorrectWithInsert(r1, r2, bestInsert);
					r2.reverseComplement();
				}
			}else if(bestInsert==RET_AMBIG){ambiguousCount++;}
			else if(bestInsert==RET_SHORT){tooShortCount++;}
			else if(bestInsert==RET_LONG){tooLongCount++;}
			else if(bestInsert==RET_NO_SOLUTION){noSolutionCount++;}
			
			r1.setInsert(bestInsert);
			
			if(findAdapterSequence && bestInsert>0){
				storeAdapterSequence(r1, bestInsert);
				r2.reverseComplement();
				storeAdapterSequence(r2, bestInsert);
				r2.reverseComplement();
			}
			
			if(originals!=null && (!ecco || bestInsert<1)){
				r1.bases=originals[0];
				r1.quality=originals[1];
				r2.bases=originals[2];
				r2.quality=originals[3];
			}
			
			if(trimReadsByOverlap && bestInsert>0){
				int trimLim=bestInsert-1;
				if(trimLim<r1.length()){
					if(verbose){System.err.println("Overlap right trimming r1 to "+0+", "+(trimLim));}
					int x=TrimRead.trimToPosition(r1, 0, trimLim, 1);
					if(verbose){System.err.println("Trimmed "+x+" bases: "+new String(r1.bases));}
				}
				if(trimLim<r2.length()){
					if(verbose){System.err.println("Overlap right trimming r2 to "+0+", "+(trimLim));}
					int x=TrimRead.trimToPosition(r2, 0, trimLim, 1);
					if(verbose){System.err.println("Trimmed "+x+" bases: "+new String(r2.bases));}
				}
			}
			
			if(bestInsert>0 || MIX_BAD_AND_GOOD){
				if(listg!=null){
					if(joined!=null){
						listg.add(joined);
					}else{
						listg.add(r1);
					}
				}
			}else if(listb!=null){
				listb.add(r1);
			}
			return bestInsert;
		}
		
		private final int preprocess(final Read r1, final Read r2, boolean qtrim){
			assert(r1!=null);
			if(!r1.validated()){r1.validate(true);}
			if(r2==null){return RET_BAD;}
			if(!r2.validated()){r2.validate(true);}
			
			if(iupacToN){
				if(r1!=null){r1.convertUndefinedTo((byte)'N');}
				if(r2!=null){r2.convertUndefinedTo((byte)'N');}
			}
			
			if(recalibrateQuality){
				CalcTrueQuality.recalibrate(r1);
				CalcTrueQuality.recalibrate(r2);
			}
			
			pairsProcessed++;
			
			if(forceTrimLeft>0 || forceTrimRight>0 || forceTrimModulo>0 || forceTrimRight2>0){
				if(r1!=null && !r1.discarded()){
					final int len=r1.length();
					final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
					final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
					final int b1=forceTrimRight>0 ? forceTrimRight : len;
					final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
					final int b=Tools.min(b0, b1, b2);
					final int x=TrimRead.trimToPosition(r1, a, b, 1);
				}
				if(r2!=null && !r2.discarded()){
					final int len=r2.length();
					final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
					final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
					final int b1=forceTrimRight>0 ? forceTrimRight : len;
					final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
					final int b=Tools.min(b0, b1, b2);
					final int x=TrimRead.trimToPosition(r2, a, b, 1);
				}
			}
			
			if(qtrim){qtrim(r1, r2, 0);}
			
			if(tadpole!=null && extendRight1>0){
				extendAndMerge(r1, r2, extendRight1, 1, false);
			}

			final int len1=r1.length(), len2=r2.length();
			
			if(len1<minReadLength && len2<minReadLength){
				return RET_BAD;
			}else if(len1<2 || len2<2){
				return RET_AMBIG;
			}
			
			if(r1.quality!=null || r2.quality!=null){
				if(minAvgQuality>0){
					if(r1.avgQuality(false, minAvgQualityBases)<minAvgQuality || r2.avgQuality(false, minAvgQualityBases)<minAvgQuality){
						//Failed quality filter
						return RET_BAD;
					}
				}
				if(useMEEfilter && useQuality){
					int maxBasesToConsider=Tools.min(Tools.max(len1, len2), len1+len2-minInsert);
					if(r1.expectedTipErrors(false, maxBasesToConsider)>maxExpectedErrors || r2.expectedTipErrors(false, maxBasesToConsider)>maxExpectedErrors){
						//Failed MEEFilter
						return RET_BAD;
					}
				}
			}
			return 1;
		}
		
		private final int extendAndMerge(Read r1, Read r2, int amt, int iters, boolean merge){
			assert(iters>0);
			assert(merge || iters==1);
			
			
			int sum1=0, sum2=0, attempted=0;
			int bestInsert=RET_AMBIG;
			for(int i=0; i<iters && (bestInsert==RET_AMBIG || bestInsert==RET_NO_SOLUTION); i++){
				
				int e1=(sum1==attempted ? extendRead(r1, amt) : 0);
				r2.reverseComplement();
				int e2=(sum2==attempted ? extendRead(r2, amt) : 0);
				r2.reverseComplement();
				
				attempted+=amt;
				sum1+=e1;
				sum2+=e2;
				
				if(merge){
					if(e1>0 || e2>0){
						bestInsert=processReadPair_inner(r1, r2);
					}else{
						break;
					}
				}
			}
			//Todo: un-extend.
			
			extensionsAttemptedT+=2;
			
			if(sum1==attempted){fullyExtendedT++;}
			else if(sum1>0){partlyExtendedT++;}
			else{notExtendedT++;}
			
			if(sum2==attempted){fullyExtendedT++;}
			else if(sum2>0){partlyExtendedT++;}
			else{notExtendedT++;}
			
			return bestInsert;
		}
		
		private final int extendRead(Read r, int amt){
			bb.clear();
			bb.append(r.bases);
			final int initialLen=r.length();
			final int extension=tadpole.extendToRight2(bb, leftCounts, rightCounts, amt, false);
			
//			extensionsAttemptedT++;
//			if(extension==amt){
//				fullyExtendedT++;
//			}else if(extension>0){
//				partlyExtendedT++;
//			}else{
//				notExtendedT++;
//			}
			
			if(extension>0){
				r.bases=bb.toBytes();
				if(r.quality!=null){
					r.quality=Arrays.copyOf(r.quality, r.bases.length);
					for(int i=initialLen; i<r.quality.length; i++){
						r.quality[i]=qfake;
					}
				}
			}
			return extension;
		}
		
		private final int mateByOverlap_ratioMode(Read r1, Read r2, int minOverlap){
			assert(useRatioMode);
			int min0=MIN_OVERLAPPING_BASES_0-MIN_OVERLAPPING_BASES_RATIO_REDUCTION;
			int min=minOverlap-MIN_OVERLAPPING_BASES_RATIO_REDUCTION;
			int x=-1;
			rvector[4]=0;

			float ratioMargin=RATIO_MARGIN;
			float maxRatio=MAX_RATIO;

			boolean overlapped=false;
			if(overlapUsingQuality && r1.quality!=null && r2.quality!=null){
				overlapped=true;
				x=BBMergeOverlapper.mateByOverlapRatio(r1, r2, aprob, bprob, rvector, 
						min0, min, minInsert0, minInsert, maxRatio, ratioMargin, RATIO_OFFSET, 0.95f, 0.95f, true);
				if(verbose){System.err.println("Result from ratiomode1:  \tinsert="+x+", bad="+rvector[2]+", ambig="+(rvector[4]==1));}
			}
			if(!overlapped || (overlapWithoutQuality && (x<0 || rvector[4]==1))){
				x=BBMergeOverlapper.mateByOverlapRatio(r1, r2, aprob, bprob, rvector, min0, min, 
						minInsert0, minInsert, maxRatio, ratioMargin, RATIO_OFFSET, 0.95f, 0.95f, false);
				if(verbose){System.err.println("Result from ratiomode2:  \tinsert="+x+", bad="+rvector[2]+", ambig="+(rvector[4]==1));}
			}
			return x;
		}
		
		private final int mateByOverlap_normalMode(Read r1, Read r2, int minOverlap){
			final int len1=r1.length(), len2=r2.length();
			boolean ambigNM=false;
			int bestInsertNM=-1;
			int bestBadNM=999999;
			
			assert(QUAL_ITERS>0);
			final int maxQualIters=(r1.quality==null || r2.quality==null ? 1 : QUAL_ITERS);
			final int maxTrims=(r1.quality==null || r2.quality==null ? 0 : TRIM_ON_OVERLAP_FAILURE);

			for(int i=0; i<maxQualIters && bestInsertNM<0 /*&& !ambigNM*/; i++){
				
				int x=BBMergeOverlapper.mateByOverlap(r1, r2, aprob, bprob, rvector, MIN_OVERLAPPING_BASES_0-i, minOverlap+i, 
						minInsert0, MISMATCH_MARGIN, MAX_MISMATCHES0, MAX_MISMATCHES, (byte)(MIN_QUALITY-2*i));
				if(x>-1){
					bestInsertNM=x;
					bestBadNM=rvector[2];
					ambigNM=(rvector[4]==1);
					break;
				}
			}


			if(loose && bestInsertNM<0){//TODO check for estimated number of overlap errors
				int x=BBMergeOverlapper.mateByOverlap(r1, r2, aprob, bprob, rvector, MIN_OVERLAPPING_BASES_0-1, minOverlap+2, 
						minInsert0, MISMATCH_MARGIN, MAX_MISMATCHES0+1, MAX_MISMATCHES+1, MIN_QUALITY-1);
				if(x>-1){
					bestInsertNM=x;
					bestBadNM=rvector[2];
					ambigNM=(rvector[4]==1);
				}
			}

			for(int trims=0, q=trimq[0]; trims<maxTrims && !qtrim1 && bestInsertNM<0 /*&& !ambigNM*/; trims++, q+=8){
				int tr1=TrimRead.trimFast(r1, false, true, q, 1+len1*4/10); //r1.length());
				int tr2=TrimRead.trimFast(r2, true, false, q, 1+len2*4/10); //r2.length());
				if(tr1>0 || tr2>0){
					int x=BBMergeOverlapper.mateByOverlap(r1, r2, aprob, bprob, rvector, MIN_OVERLAPPING_BASES_0-1, minOverlap, 
							minInsert0, MISMATCH_MARGIN, MAX_MISMATCHES0, MAX_MISMATCHES, MIN_QUALITY);
					if(x>-1){
						bestInsertNM=x;
						bestBadNM=rvector[2];
						ambigNM=(rvector[4]==1);
						trims=maxTrims;
					}
				}
			}
			if(verbose){System.err.println("Result from normalmode:  \tinsert="+bestInsertNM+", bad="+bestBadNM+", ambig="+ambigNM);}
			
			rvector[0]=bestInsertNM;
			rvector[2]=bestBadNM;
			rvector[4]=(ambigNM ? 1 : 0);
			return bestInsertNM;
		}
		
		private final int calcMinOverlapFromEntropy(final Read r1, final Read r2){
			if(!useEntropy){return MIN_OVERLAPPING_BASES;}
			final int minOverlap;
			if(loose){
				final int len1=r1.length(), len2=r2.length();
				int a=BBMergeOverlapper.calcMinOverlapByEntropy(r1.bases, entropyK, kmerCounts, minEntropyScore);
				int b=BBMergeOverlapper.calcMinOverlapByEntropy(r2.bases, entropyK, kmerCounts, minEntropyScore);
				float errorRate=r1.expectedErrors(false, len1)/len1+r2.expectedErrors(false, len2)/len2;
				minOverlap=(int)(Tools.max(MIN_OVERLAPPING_BASES, Tools.max(a, b))*(1+Tools.min(0.04f, errorRate)*4f));
			}else{
				int a=BBMergeOverlapper.calcMinOverlapByEntropyTail(r1.bases, entropyK, kmerCounts, minEntropyScore);
				int b=BBMergeOverlapper.calcMinOverlapByEntropyHead(r2.bases, entropyK, kmerCounts, minEntropyScore);
				minOverlap=Tools.max(MIN_OVERLAPPING_BASES, Tools.max(a, b));
			}
			return minOverlap;
		}
		
		private final int lookForAdapters(final Read r1, final Read r2){
			assert(lowercaseAdapters);
			if(!lowercaseAdapters){return -1;}
			if(!Character.isLowerCase(r1.bases[r1.length()-1]) || !Character.isLowerCase(r2.bases[0])){return -1;}
			
			final int lower1=r1.trailingLowerCase(), lower2=r2.leadingLowerCase();

			final int upper1=r1.length()-lower1, upper2=r2.length()-lower2;
			final int newlen=Tools.min(upper1, upper2);
			int good=0, bad=0;

			for(int i=0; i<newlen; i++){
				byte a=r1.bases[i];
				byte b=r2.bases[i+lower2];
				if(a!='N' && b!='N'){
					if(a==b){good++;}
					else{bad++;}
				}
			}
			if(bad*4<=good){
				rvector[0]=newlen;
				rvector[2]=bad;
				rvector[4]=0;
				return newlen;
			}
			return -1;
		}
		
		private final int mateByOverlap(Read r1, Read r2){
			final int len1=r1.length(), len2=r2.length();
			
			final int minOverlap=calcMinOverlapFromEntropy(r1, r2);
			if(verbose){System.err.println("minOverlap: "+minOverlap);}
			
			//TODO: Currently this is not used for anything.
			final int bestInsertAD;
			final int bestBadAD;
			if(lowercaseAdapters){
				bestInsertAD=lookForAdapters(r2, r2);
				bestBadAD=(bestInsertAD>=0 ? rvector[2] : 0);
			}
			
			if(aprob==null || aprob.length<Tools.max(len1, len2)){aprob=new float[Tools.max(len1, len2)];}
			if(bprob==null || bprob.length<Tools.max(len1, len2)){bprob=new float[Tools.max(len1, len2)];}
			
			final boolean ambigRM;
			final int bestBadRM, bestInsertRM;
			if(useRatioMode){
				bestInsertRM=mateByOverlap_ratioMode(r1, r2, minOverlap);
				bestBadRM=rvector[2];
				ambigRM=(bestInsertRM>-1 ? rvector[4]==1 : false);
			}else{
				bestInsertRM=-1;
				bestBadRM=0;
				ambigRM=false;
			}
			
			final boolean ambigNM;
			final int bestInsertNM, bestBadNM;
			if(useNormalMode && ((!requireRatioMatch && (bestInsertRM<0 || ambigRM)) || (requireRatioMatch && (bestInsertRM>0 && !ambigRM)))){
				bestInsertNM=mateByOverlap_normalMode(r1, r2, minOverlap);
				bestBadNM=rvector[2];
				ambigNM=(bestInsertNM>-1 ? rvector[4]==1 : false);
			}else{
				ambigNM=false;
				bestInsertNM=-1;
				bestBadNM=99999;
			}
			
			boolean ambig;
			int bestBad, bestInsert;
			if(requireRatioMatch && useNormalMode && useRatioMode){
				ambig=ambigRM || ambigNM;
				bestBad=bestBadRM;
				bestInsert=(bestInsertNM==bestInsertRM ? bestInsertNM : -1);

				if(verbose){System.err.println("Result after rrm:  \tinsert="+bestInsertNM+", bad="+bestBadNM+", ambig="+ambigNM);}
			}else if(useRatioMode && bestInsertRM>-1 && !ambigRM){
				ambig=ambigRM;
				bestBad=bestBadRM;
				bestInsert=bestInsertRM;
			}else{
				ambig=ambigNM;
				bestBad=bestBadNM;
				bestInsert=bestInsertNM;
			}
			
			if(bestBad>MAX_MISMATCHES_R){ambig=true;}
			
			if(ambig){return RET_AMBIG;}
			else if(bestInsert<0){return RET_NO_SOLUTION;}

			//TODO:  Crucial!  This block can vastly reduce merge rate, particularly if quality values are inaccurate.
			if(useQuality && r1.quality!=null && r2.quality!=null){
				if(useEfilter && bestInsert>0 && !ambig){
					float bestExpected=BBMergeOverlapper.expectedMismatches(r1, r2, bestInsert);
					if((bestExpected+efilterOffset)*efilterRatio<bestBad){ambig=true;}
					if(verbose){System.err.println("Result after efilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
				}

				if(pfilterRatio>0 && bestInsert>0 && !ambig){
					float probability=BBMergeOverlapper.probability(r1, r2, bestInsert);
					if(probability<pfilterRatio){bestInsert=-1;}
					if(verbose){System.err.println("Result after pfilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
				}
			}
			
			if(ambig){return RET_AMBIG;}
			r1.errors=bestBad;
			return bestInsert>0 ? bestInsert : RET_NO_SOLUTION;
		}
		
		/**
		 * 
		 * @param r1 Read1
		 * @param r2 Read2
		 * @return A return code (RET_)
		 */
		private final int processReadPair(final Read r1, final Read r2){
			
			{
				final int x=preprocess(r1, r2, (qtrim1 && !qtrim2));
				if(x<0){return x;}
			}
			
			if(parseCustom){
				int trueSize=-1;
				if(r1.id.startsWith("insert=")){
					trueSize=GradeMergedReads.parseInsert(r1.id);
				}else{
					r1.setMapped(true);
					r2.setMapped(true);
					trueSize=Read.insertSizeMapped(r1, r2, ignoreMappingStrand);
				}
				if(verbose){System.err.println("True Insert: "+trueSize);}
				r1.setInsert(trueSize);
			}
			
			r2.reverseComplement();
			
			byte[] qual1=r1.quality, qual2=r2.quality;
			if(!useQuality){//strip qualities
				r1.quality=r2.quality=null;
			}
			
			int bestInsert=processReadPair_inner(r1, r2);
			if(qtrim2){
				for(int iter=0; iter<trimq.length && bestInsert<0; iter++){
					r1.quality=qual1;
					r2.quality=qual2;
					//				r2.reverseComplement();
					//				qtrim(r1, r2);
					//				r2.reverseComplement();

					TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq[iter], 1);
					TrimRead.trimFast(r2, qtrimRight, qtrimLeft, trimq[iter], 1);//Reversed because read is rcomped

					qual1=r1.quality;
					qual2=r2.quality;
					bestInsert=processReadPair_inner(r1, r2);
				}
			}
			
			if(tadpole!=null){
				if(eccTadpole && (bestInsert==RET_AMBIG || bestInsert==RET_NO_SOLUTION)){
					int c1=tadpole.errorCorrect(r1);
					int c2=tadpole.errorCorrect(r2);
					if(c1>0 || c2>0){
						bestInsert=processReadPair_inner(r1, r2);
					}
				}

				if(extendRight2>0 && (bestInsert==RET_AMBIG || bestInsert==RET_NO_SOLUTION)){
					bestInsert=extendAndMerge(r1, r2, extendRight2, extendIterations, true);
				}
			}
			
			if(useKFilter && bestInsert>kmerLength){
				Read joined=r1.joinRead(bestInsert);
				if(useKFilter){
					int cov=BBMergeOverlapper.minCoverage(joined, tadpole, kmerLength, filterCutoff);
					if(cov<filterCutoff){bestInsert=RET_NO_SOLUTION;}
					if(verbose){System.err.println("Result after kfilter:  \tinsert="+bestInsert);}
				}
			}
			
			if(!useQuality){//restore qualities
				r1.quality=qual1;
				r2.quality=qual2;
			}
			r2.reverseComplement();
			return bestInsert;
		}
		
		/**
		 * 
		 * @param r1 Read1
		 * @param r2 Read2
		 * @return A return code (RET_)
		 */
		private final int processReadPair_inner(final Read r1, final Read r2){
			int bestInsert=-1;
			boolean ambig;
			
			if(USE_MAPPING && r1.chrom==r2.chrom && r1.start<r1.stop && ((r1.mapped() || r1.synthetic()) && (r2.mapped() || r2.synthetic()))){
				bestInsert=r1.insert();
				ambig=false;
			}else{
				if(MATE_BY_OVERLAP){
					bestInsert=mateByOverlap(r1, r2);
					ambig=(bestInsert==RET_AMBIG);
				}else{
					ambig=false;
					bestInsert=-1;
				}
			}
			
			if(ambig){return RET_AMBIG;}
			else if(bestInsert>0){
				if(bestInsert<minInsert){return RET_SHORT;}
				else if(bestInsert>maxReadLength){return RET_LONG;}
				return bestInsert;
			}
			else{return RET_NO_SOLUTION;}
		}
		
		private void storeAdapterSequence(Read r, int insert){
			LongList[] lists=adapterCountsT[r.pairnum()];
			byte[] bases=r.bases;
			for(int i=insert, j=0; i<bases.length; i++, j++){
				byte b=bases[i];
				int num=AminoAcid.baseToNumber[b];
				if(num>=0){
					lists[num].increment(j);
				}
			}
		}
		
		
		/*--------------------------------------------------------------*/
		
		final LongList[][] adapterCountsT=new LongList[2][4];
		
		final byte qfake=Shared.FAKE_QUAL;
		
		private final int[] rvector=new int[5];

		private final int[] rightCounts=new int[4];
		private final int[] leftCounts=(extendThroughLeftJunctions ? null : new int[4]);
		
		private final ByteBuilder bb=new ByteBuilder();
		
		final long[] hist=new long[histlen];
		final short[] kmerCounts;
		
		private float[] aprob, bprob;

		long pairsProcessed=0;
		long matedCount=0;
		long correctCount=0;
		long ambiguousCount=0;
		long tooShortCount=0;
		long tooLongCount=0;
		long incorrectCount=0;
		long noSolutionCount=0;
		long insertSumCorrect=0;
		long insertSumIncorrect=0;
		int insertMax=0;
		int insertMin=999999999;

		long fullyExtendedT=0;
		long partlyExtendedT=0;
		long notExtendedT=0;
		long extensionsAttemptedT=0;
		
		private final ConcurrentReadInputStream cris;
		private final ConcurrentReadOutputStream rosgood;
		private final ConcurrentReadOutputStream rosbad;
		private final ConcurrentReadOutputStream rosi;
		
		private final boolean joinReads;
		private final boolean trimReadsByOverlap;
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1;
	private String in2;
	
	private ArrayList<String> extra=new ArrayList<String>();

	private String out1=null;
	private String out2=null;
	private String outb1=null;
	private String outb2=null;
	private String outinsert=null;
	private String ihist=null;
	private String outAdapter=null;
	private String outCardinality=null;
	
	private final LogLog loglog;
	
	private long maxReads=-1;
	private boolean join=true;
	private boolean ecco=false;
	private boolean trimByOverlap=false;
	
	private float pfilterRatio=0.00002f;
	private float efilterRatio=6f;
	private float efilterOffset=0.05f;
	private boolean useEfilter=true;
	private boolean useMEEfilter=false;
	
	private boolean ordered=false;
	private boolean overlapUsingQuality=false;
	private boolean overlapWithoutQuality=true;
	private boolean useKFilter=false;
	private int filterCutoff=1;
	private int kmerLength=31;
	private boolean prealloc=false;
	private boolean prefilter=false;
	private ArrayList<String> extraFiles;

	private boolean useEntropy=true;
	private int entropyK=3;
	private int minEntropyScore=39;//30 loose;//39 normal;//44 strict;
	
	private long sampleseed=-1;
	private float samplerate=1;
	
	private boolean findAdapterSequence=false;
	
	private final LongList[][] adapterCounts=new LongList[2][4];
	
	private final Tadpole tadpole;
	private int extendRight1=0;
	private int extendRight2=0;
	private int extendIterations=1;
	private boolean eccTadpole=false;
	private boolean shave=false;
	private boolean rinse=false;
	
	private boolean extendThroughLeftJunctions=true;
	private int minCountSeed=3, minCountExtend=2;
	private float branchMult1=20;
	private float branchMult2=3;
	private float minProb=0.5f;
	private int branchLowerConst=3;
	
	/*--------------------------------------------------------------*/
	
	private static ThreadLocal<int[]> localRvector=new ThreadLocal<int[]>();
	
	static boolean errorState=false;

	private static boolean showFullArgs=true;

	/** Recalibrate quality scores using matrices */
	static boolean recalibrateQuality=false;
	static boolean useQuality=true;
	static boolean qtrimRight=false;
	static boolean qtrimLeft=false;
//	static boolean untrim=false;
	static byte[] trimq=new byte[] {6};
	static byte minAvgQuality=0;
	static int minAvgQualityBases=0;
	static float maxExpectedErrors=0;
	static int minReadLength=1;
	static int maxReadLength=-1;
	static int minInsert=35;
	static int minInsert0=-1;
	static boolean qtrim1=false;
	static boolean qtrim2=false;
	static int TRIM_ON_OVERLAP_FAILURE=1;
	static int QUAL_ITERS=3;
	static boolean parseCustom=false;
	
	static int forceTrimLeft;
	static int forceTrimRight;
	static int forceTrimRight2;
	/** Trim right bases of the read modulo this value. 
	 * e.g. forceTrimModulo=50 would trim the last 3bp from a 153bp read. */
	static int forceTrimModulo;
	
	static boolean strict=false;
	static boolean vstrict=false;
	static boolean ustrict=false;
	static boolean xstrict=false;
	static boolean loose=false;
	static boolean vloose=false;
	static boolean uloose=false;
	static boolean xloose=false;
	static boolean fast=false;
	
	/** If true, interpret lowercase bases as adapter sequence */
	static boolean lowercaseAdapters=false;
	
	private static final int histlen=2000; 
	static long[] histTotal=new long[histlen];
	static int bin=1;

	static long readsProcessedTotal=0;
	static long matedCountTotal=0;
	static long correctCountTotal=0;
	static long ambiguousCountTotal=0;
	static long tooShortCountTotal=0;
	static long tooLongCountTotal=0;
	static long incorrectCountTotal=0;
	static long noSolutionCountTotal=0;
	static long insertSumCorrectTotal=0;
	static long insertSumIncorrectTotal=0;
	static long fullyExtendedTotal=0;
	static long partlyExtendedTotal=0;
	static long notExtendedTotal=0;
	static long extensionsAttempted=0;
	
	static int insertMinTotal=999999999;
	static int insertMaxTotal=0;
	
	private static int MIN_OVERLAPPING_BASES=11;
	private static int MIN_OVERLAPPING_BASES_0=8;
	private static int MISMATCH_MARGIN=2;
	private static int MIN_OVERLAPPING_BASES_RATIO_REDUCTION=3;
	
	static boolean useRatioMode=true;
	static boolean useNormalMode=false;
	static boolean requireRatioMatch=false;
	static int MAX_MISMATCHES_R=20;
	static float MAX_RATIO=0.09f;
	static float RATIO_MARGIN=5.5f;
	static float RATIO_OFFSET=0.55f;
	
	public static int MAX_MISMATCHES=3;
	public static int MAX_MISMATCHES0=3;
	public static byte MIN_QUALITY=10;
	
	public static final int RET_NO_SOLUTION=-1;
	public static final int RET_AMBIG=-2;
	public static final int RET_BAD=-3;
	public static final int RET_SHORT=-4;
	public static final int RET_LONG=-5;
	
	/** Skip alignment and calculate insert from mapping info */ 
	protected static boolean USE_MAPPING=false;
	protected static final boolean ignoreMappingStrand=false;
	
	private static boolean MATE_BY_OVERLAP=true;
//	private static boolean SKIP_MATED_READS=false;
//	private static boolean OUTPUT_FAILED=true;
	private static boolean MIX_BAD_AND_GOOD=false;
	private static boolean NONZERO_ONLY=true;
	private static boolean showHistStats=true;
	
	private static boolean overwrite=true;
	private static boolean append=false;
	private static final boolean verbose=false;
	
	private static boolean iupacToN=false;
	
	private static int THREADS=-1;
	private static String version="8.9";
	
}
