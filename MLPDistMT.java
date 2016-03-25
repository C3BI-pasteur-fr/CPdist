import java.io.*;
import java.util.*;
import java.text.*;

// native code compilation : gcj -Ofast --no-bounds-check -march=native --main=MLPdist MLPdist.java -o mlpdist

// USAGE : MLPdist -i data.txt  -t outlierThreshold  -m [P,JC69,F81]

public class MLPDistMT {

    // constants
    final static String NO_FILE = "NOFILE";
    final static String NO_SEQ = "?";
    final static String BLANK = "                                                                                                                                                                                                                                                                                                                         ";
    final static double SMALL = 1E-10;


    // options
    static File datafile;              // -i : input file containing all FASTA-formatted MSAs
    static File outfile;               // -o : output file containing the distance matrix
    static double outlierThreshold;    // -t : default 1%
    static double precision;           // -p : default 2%
    static int minlocus;               // -n : any isolate shared by < minlocus is discarded
    static String method;              // -m : default P ; otherwise JC69, F81
    static int bootstrap;              // -b : default 0
    static int core;                   // -c : no. threads to use; default 1 
    static boolean verbose;            // -v

    // io
    static NumberFormat df;
    static BufferedReader in;
    static BufferedWriter out;

    // data
    static int k;                                  // no. MSAs
    static ArrayList<File> infile;                 // MSA file names
    static int n;                                  // no. taxa
    static ArrayList<String> tax;                  // taxon names
    static int mtl;                                // max tax length
    static boolean[][] exists;                     // existence matrix : exists[i][m]
    static short[][] seqid;                        // sequence identifier table : seqid[i][m] has the sequence seqset.get(m).get(i)
    static ArrayList<ArrayList<String>> seqset;    // sequence set : seqset.get(m) contains every distinct aligned sequence for locus m
    static Point2D[][][] pt2m;                     // (s,l) point matrix : pt2d[i][j][m] where j < i < seqset.get(m).size()
    static int[] cno;                              // char no. of each MSA

    // NB estimate
    static Point2d[][] nbprm;
    static NBfactory nbf;
    static ArrayList<String> keya; 
    static ArrayList<Point2D> loca;
       
    // multi-threading
    static NBfactory[] nbfa;
    static Thread[] thread;
    
    // pairwise distances
    static double[][] dm;

    // bootstrap
    static Random rand;
    static int[] bid;    

    // stuffs
    static int m, i, j, c, s, l, o, t, it, jt, b, mb;
    static double up, down, p, r, pv;
    static String line, lbl, si, sj, lv;
    static char ci, cj;
    static short im, jm;
    static boolean ok;
    static File f;
    static StringBuffer sb;
    static StringBuilder sB;
    static ArrayList<Integer> inta;
    static Point2D sl;


    public static void main(String[] args) throws IOException {



	// ######################################
	// ##### init parameters            #####
	// ######################################
	df = NumberFormat.getNumberInstance(Locale.US); df.setGroupingUsed(false); df.setMaximumFractionDigits(8); df.setMinimumFractionDigits(8);
	datafile = new File(NO_FILE);
	outfile = new File(NO_FILE);
	outlierThreshold = 0.01;
	precision = 1E-4;
	method = "P";
	bootstrap = 0;
	minlocus = 0; 
	verbose = false;
	core = 1;


	// ######################################
	// ##### reading options            #####
	// ######################################
	o = -1;
	while ( ++o < args.length ) {
	    if ( args[o].equals("-i") ) {
		datafile = new File(args[++o]);
		if ( ! datafile.exists() ) { System.out.println(args[0] + " does not exists (option -i)"); System.exit(0); }
		continue; 
	    }
	    if ( args[o].equals("-o") ) { outfile = new File(args[++o]);; continue; }
	    if ( args[o].equals("-t") ) { outlierThreshold = Double.parseDouble(args[++o]); continue; }
	    if ( args[o].equals("-m") ) { method = args[++o]; continue; }
	    if ( args[o].equals("-n") ) { minlocus = Integer.parseInt(args[++o]); continue; }
	    if ( args[o].equals("-b") ) { bootstrap = Integer.parseInt(args[++o]); continue; }
	    if ( args[o].equals("-c") ) { core = Integer.parseInt(args[++o]); continue; }
	    if ( args[o].equals("-v") ) { verbose = true; continue; }
	}
	if ( datafile.toString().equals(NO_FILE) ) { System.out.println(" no input file (option -i)"); System.exit(0); }
	if ( outfile.toString().equals(NO_FILE) ) outfile = new File(datafile.toString() + ".d");


	// ######################################
	// ##### reading datafile           #####
	// ######################################
	infile = new ArrayList<File>();
	in = new BufferedReader(new FileReader(datafile));
	while ( true ) {
	    try { if ( (line=in.readLine().trim()).length() == 0 ) continue; } catch ( NullPointerException e ) { in.close(); break; }
	    if ( ! (f=new File(line)).exists() ) { System.out.println("WARNING: " + line + " does not exist"); continue; }
	    infile.add(f);
	}
	k = infile.size();
	System.err.println(k + " input files    ");


	// #################################################
	// ##### reading infiles to get taxon names    #####
	// #################################################
	tax = new ArrayList<String>(); inta = new ArrayList<Integer>(); mtl = 0; m = -1;
	while ( ++m < k ) {
	    in = new BufferedReader(new FileReader(infile.get(m)));
	    while ( true ) {
		try { if ( ! (line=in.readLine().trim()).startsWith(">") ) continue; } catch ( NullPointerException e ) { in.close(); break; }
		lbl = line.substring(1).trim();
		if ( (i=tax.indexOf(lbl)) < 0 ) { tax.add(lbl); inta.add(new Integer(1)); mtl = ( mtl < (l=lbl.length()) ) ? l : mtl; }
		else inta.set(i, new Integer(inta.get(i).intValue()+1));
		
	    }
	}
	n = tax.size(); // System.out.println(tax.toString());
	if ( minlocus > 0 ) { i = n; while ( --i >= 0 ) if ( inta.get(i).intValue() < minlocus ) { tax.remove(i); inta.remove(i); } }
	n = tax.size(); // System.out.println(tax.toString());
	System.err.println(n + " taxa");


	// ##########################################################
	// ##### reading infiles and storing aligned sequence   #####
	// ##########################################################
	seqid = new short[n][k]; exists = new boolean[n][k]; seqset = new ArrayList<ArrayList<String>>(k); pt2m = new Point2D[k][][]; cno = new int[k]; i = n; while ( --i >= 0 ) Arrays.fill(exists[i], false);
	m = -1;
	while ( ++m < k ) {
	    lv = "[" + m + "/" + k + "] processing files..."; if ( verbose ) System.out.print(lv); else if ( m % 10 == 0 ) System.err.print("\r" + BLANK.substring(0, lv.length()) + "\r" + lv); 
	    // ##### reading fasta file #####
	    i = -1; sb = new StringBuffer(""); seqset.add( new ArrayList<String>() ); cno[m] = Integer.MAX_VALUE; in = new BufferedReader(new FileReader(infile.get(m))); 
	    while ( true ) {
		try { line = in.readLine().trim(); } catch ( NullPointerException e ) { in.close(); break; }
		if ( line.startsWith(">") ) {
		    if ( i != -1 ) { if ( (j=seqset.get(m).indexOf(si=sb.toString().toUpperCase())) == -1 ) { seqset.get(m).add(si); cno[m] = Math.min(cno[m], si.length()); j = seqset.get(m).size(); --j; } seqid[i][m] = (short) j; exists[i][m] = true; }
		    i = tax.indexOf(line.substring(1).trim()); sb = new StringBuffer(sb.length());
		    continue;
		}
		sb = sb.append(line);
	    }
	    if ( i != -1 ) { if ( (j=seqset.get(m).indexOf(si=sb.toString().toUpperCase())) == -1 ) { seqset.get(m).add(si); cno[m] = Math.min(cno[m], si.length()); j = seqset.get(m).size(); --j; } seqid[i][m] = (short) j; exists[i][m] = true; } 
	    //i = -1; while ( ++i < n ) if ( exists[i][m] ) System.out.println(i + " " + (tax.get(i) + BLANK).substring(0,mtl) + " " + (si=seqa[i][m]).substring(0, Math.min(100, si.length()))); else System.out.println((tax.get(i) + BLANK).substring(0,mtl) + " -");
	    // ##### estimating locus pairwise distances #####
	    pt2m[m] = new Point2D[i=seqset.get(m).size()][];
	    while ( --i >= 0 ) {
		pt2m[m][i] = new Point2D[++i]; 
		c = (l = (si=seqset.get(m).get(--i)).length()); while ( --c >= 0 ) switch ( (ci=si.charAt(c)) ) { case '-': case '?': --l; continue; } pt2m[m][i][i] = new Point2D(0, l); 
		j = i; 
		while ( --j >= 0 ) {
		    c = cno[m]; sj = seqset.get(m).get(j);
		    s = 0; l = c; while ( --c >= 0 ) switch ( (ci=si.charAt(c)) ) { case '-': case '?': --l; continue; default: switch ( (cj=sj.charAt(c)) ) { case '-': case '?': --l; continue; default: s += ( ci != cj ) ? 1 : 0; continue; } } 
		    pt2m[m][i][j] = new Point2D(s, l);
		}
	    }
	    if ( verbose ) System.out.print("\r" + BLANK.substring(0, lv.length()) + "\r"); 
	}
	lv = "[" + m + "/" + k + "] processing files..."; if ( ! verbose ) System.err.print("\r" + BLANK.substring(0, lv.length()) + "  \r");



	// ###########################################################################################
	// ##### estimating NB parameters and estimating pairwise distances                      #####
	// ###########################################################################################
	keya = new ArrayList<String>(n*(n-1)/2); loca = new ArrayList<Point2D>(n*(n-1)/2); sB = new StringBuilder(2*k); //m = k; while ( --m >= 0 ) sB = sB.append(' ');
	nbprm = new Point2d[n][]; dm = new double[n][]; i = n; while ( --i >= 0 ) { nbprm[i] = new Point2d[i]; dm[i] = new double[i]; }
	o = n*(n-1)/2; c = 0;
	if ( (core == 1) || (o <= core) ) { // one-thread computing
	    i = n;
	    while ( --i >= 0 ) {
		j = i;
		while ( --j >= 0 ) {
		    //++c; nbf = new NBfactory(new Point2D(i, j), seqid, exists, pt2m, outlierThreshold); nbf.run(); nbprm[i][j] = nbf.getNBparameters();
		    
		    //++c; sb = new StringBuffer(2*k); m = k; while ( --m >= 0 ) sb = sb.append( (exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0 )) ? sl.getX() : "?" ).append(" ");
		    ++c; sB = sB.delete(0, 2*k); m = k; while ( --m >= 0 ) { ok = exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0); sB = sB.append((char) ((ok) ? 33+sl.getX() : 32)).append((char) ((ok) ? cno[m]+33-sl.getY() : 32)); }
		    b = Collections.binarySearch(keya, (line=sB.toString())); //System.out.println(b + " " + keya.size());
		    if ( b < 0 ) { b = -(++b); keya.add(b, line); loca.add(b, new Point2D(i, j) ); nbf = new NBfactory(new Point2D(i, j), seqid, exists, pt2m, outlierThreshold); nbf.run(); nbprm[i][j] = nbf.getNBparameters(); }
		    else nbprm[i][j] = nbprm[loca.get(b).getX()][loca.get(b).getY()].copy();
		    
		    p = nbprm[i][j].getX(); r = nbprm[i][j].getY(); up = 0; down = 0; m = k; 
		    while ( --m >= 0 ) if ( exists[i][m] && exists[j][m] && ( ( sl = ((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm].copy() : pt2m[m][jm][im].copy() ).getY() > 0 ) ) { pv = Math.min( (pv=nbcdf((s=sl.getX()), p, r*(l=sl.getY()))) , 1-pv+nbpmf(s, p, r*l) ); up += pv*s; down += pv*l; }
		    dm[i][j] = up / down;

		    if ( verbose ) {
			up = 0; down = 0; m = k;
			while ( --m >= 0 )
			    if ( exists[i][m] && exists[j][m] && ( ( sl = ( (im=seqid[i][m]) > (jm=seqid[j][m]) ) ? pt2m[m][im][jm].copy() : pt2m[m][jm][im].copy() ).getY() > 0 ) ) {
				up += sl.getX(); down += sl.getY();
			    }
			System.out.println("["+c+"/"+o+"] ("+tax.get(i)+" "+tax.get(j)+") ("+df.format(p)+" "+df.format(r)+") " + df.format(up/down) + " " + df.format(p*r/(1-p)) + " " + df.format(dm[i][j])); 
		    }
		    else if ( c % 50 == 0 ) System.err.print("\r                                           \r[" + (100*c/o) + "%] estimating pairwise distances ...");
		}
	    }
	}
	else { // multi-thread computing
	    nbfa = new NBfactory[core]; thread = new Thread[core]; 
	    i = 0; j = -1; t = core; 
	    while ( --t >= 0 ) {
		j = ( (++j) == i ) ? 0 : j; i += ( j == 0 ) ? 1 : 0; 
		//++c; dm[i][j] = -1; nbfa[t] = new NBfactory(new Point2D(i, j), seqid, exists, pt2m, outlierThreshold); thread[t] = new Thread(nbfa[t]); thread[t].start();  
		sB = sB.delete(0, 2*k); m = k; while ( --m >= 0 ) { ok = exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0); sB = sB.append((char) ((ok) ? 33+sl.getX() : 32)).append((char) ((ok) ? cno[m]+33-sl.getY() : 32)); }
		++c; dm[i][j] = -1; 
		//sb = new StringBuffer(2*k); m = k; while ( --m >= 0 ) sb = sb.append( (exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0 )) ? sl.getX() : "?" ).append(" ");
		//m = k; while ( --m >= 0 ) sB.setCharAt(m , (char) ((exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0 )) ? 33+sl.getX() : 32) );
		
		b = Collections.binarySearch(keya, (line=sB.toString())); 
		if ( b < 0 ) { b = -(++b); keya.add(b, line); loca.add(b, new Point2D(i, j) ); nbfa[t] = new NBfactory(new Point2D(i, j), seqid, exists, pt2m, outlierThreshold); thread[t] = new Thread(nbfa[t]); thread[t].start(); }
		else ++t;

	    }
	    do {
		t = core;
		while ( --t >= 0 ) {
		    if ( thread[t] == null ) { // launching new thread
			j = ( (++j) == i ) ? 0 : j; i += ( j == 0 ) ? 1 : 0; 
			//if ( i < n ) { ++c; dm[i][j] = -1; nbfa[t] = new NBfactory(new Point2D(i, j), seqid, exists, pt2m, outlierThreshold); thread[t] = new Thread(nbfa[t]); thread[t].start(); }

			if ( i < n ) {
			    ++c; dm[i][j] = -1; 
			    //sb = new StringBuffer(2*k); m = k; while ( --m >= 0 ) sb = sb.append( (exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0 )) ? sl.getX() : "?" ).append(" ");
			    //m = k; while ( --m >= 0 ) sB.setCharAt(m, (char) ((exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0 )) ? 33+sl.getX() : 32) );
			    //m = k; while ( --m >= 0 ) sB = sB.append((char) ((ok = (exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0))) ? 33+sl.getX() : 32)).append((char) (( ok ) ? cno[m]+33-sl.getY() : 32));
			    sB = sB.delete(0, 2*k); m = k; while ( --m >= 0 ) { ok = exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0); sB = sB.append((char) ((ok) ? 33+sl.getX() : 32)).append((char) ((ok) ? cno[m]+33-sl.getY() : 32)); }
			    
			    b = Collections.binarySearch(keya, (line=sB.toString())); //System.out.println(line);
			    if ( b < 0 ) { b = -(++b); keya.add(b, line); loca.add(b, new Point2D(i, j) ); nbfa[t] = new NBfactory(new Point2D(i, j), seqid, exists, pt2m, outlierThreshold); thread[t] = new Thread(nbfa[t]); thread[t].start(); }
			    else { ++t; if ( verbose ) System.out.println("["+c+"/"+o+"] ("+i+","+j+") ("+tax.get(i)+" "+tax.get(j)+") already being estimated"); }
			}
			
			continue;
		    }

		    if ( thread[t].getState().equals(Thread.State.TERMINATED) ) {  // getting NBfactory results
			// getting taxon pair and NB parameters
			nbprm[it=nbfa[t].getIndexes().getX()][jt=nbfa[t].getIndexes().getY()] = nbfa[t].getNBparameters(); p = nbprm[it][jt].getX(); r = nbprm[it][jt].getY();
			nbfa[t] = null; thread[t] = null; 
			// computing pairwise distance
			up = 0; down = 0; m = k; 
			while ( --m >= 0 ) if ( exists[it][m] && exists[jt][m] && ((sl = ((im=seqid[it][m])>(jm=seqid[jt][m])) ? pt2m[m][im][jm].copy() : pt2m[m][jm][im].copy() ).getY() > 0) ) { pv = Math.min((pv=nbcdf((s=sl.getX()), p, r*(l=sl.getY()))), 1-pv+nbpmf(s, p, r*l) ); up += pv*s; down += pv*l; }
			dm[it][jt] = up / down;
			if ( verbose ) {
			    up = 0; down = 0; m = k; while ( --m >= 0 ) if ( exists[it][m] && exists[jt][m] && ( ( sl = ( (im=seqid[it][m]) > (jm=seqid[jt][m]) ) ? pt2m[m][im][jm].copy() : pt2m[m][jm][im].copy() ).getY() > 0 ) ) { up += sl.getX(); down += sl.getY(); }
			    System.out.println("["+c+"/"+o+"] ("+it+","+jt+") ("+tax.get(it)+" "+tax.get(jt)+") ("+df.format(p)+" "+df.format(r)+") " + df.format(up/down) + " " + df.format(p*r/(1-p)) + " " + df.format(dm[it][jt])); 
			}
			else if ( c % 50 == 0 ) System.err.print("\r                                           \r[" + (100*c/o) + "%] estimating pairwise distances ...");
		    }
		}
		while ( Thread.activeCount() > core ) { }                                     // waiting for some threads to terminate
		if ( c == o ) {                                                               // all pairs were considered
		    b = 0; t = core; while ( --t >= 0 ) b += ( thread[t] == null ) ? 1 : 0;   // if b < core, then some threads are already running
		    if ( b < core ) { t = 0; while ( t >= 0 ) { t = core; while ( --t >= 0 ) if ( (thread[t] != null) && ! thread[t].getState().equals(Thread.State.TERMINATED) ) break; } } else ++c;
		}
	    } while ( c <= o );

	    // filling the empty entries of the distance matrix, i.e. here each dm[i][j] == -1 corresponds to NB parameters that were already estimated
	    i = n;
	    while ( --i >= 0 ) {
		System.err.print("\r                                           \r[" + (100*(n-i)/n) + "%] finishing ...");
		j = i;
		while ( --j >= 0 ) {
		    if ( dm[i][j] != -1 ) continue;
		    //sb = new StringBuffer(2*k); m = k; while ( --m >= 0 ) sb = sb.append( (exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0 )) ? sl.getX() : "?" ).append(" ");
		    //m = k; while ( --m >= 0 ) sB.setCharAt(m, (char) ((exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0 )) ? 33+sl.getX() : 32) );
		    sB = sB.delete(0, 2*k); m = k; while ( --m >= 0 ) { ok = exists[i][m] && exists[j][m] && ((sl=((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im]).getY() > 0); sB = sB.append((char) ((ok) ? 33+sl.getX() : 32)).append((char) ((ok) ? cno[m]+33-sl.getY() : 32)); }
		    
		    b = Collections.binarySearch(keya, (line=sB.toString())); 
		    nbprm[i][j] = nbprm[it=loca.get(b).getX()][jt=loca.get(b).getY()].copy();
		    dm[i][j] = dm[it][jt];
		}
	    }
	}
    

    

	// ###########################################################################################
	// ##### outputing distance matrix                                                       #####
	// ###########################################################################################
	df.setMaximumFractionDigits(12); df.setMinimumFractionDigits(12);
	out = new BufferedWriter(new FileWriter(outfile));
	out.write(" " + n); out.newLine(); i = -1; 
	while ( ++i < n ) {
	    out.write((tax.get(i) + BLANK).substring(0,mtl));
	    j = -1; while ( ++j < i ) out.write(" " + df.format(dm[i][j])); out.write(" " + df.format(0)); while ( ++j < n ) out.write(" " + df.format(dm[j][i]));
	    out.newLine(); 
	}
	out.newLine(); 

	
	// ###########################################################################################
	// ##### bootstraping                                                                    #####
	// ###########################################################################################
	if ( bootstrap > 0 ) {
	    rand = new Random(bootstrap * n); bid = new int[k]; o = -1;
	    while ( ++o < bootstrap ) {
		if ( o % 5 == 0 ) System.err.print("\r                                                    \r[" + (100*o/bootstrap) + "%] estimating bootstrap pairwise distances ...");
		m = k; while ( --m >= 0 ) bid[m] = rand.nextInt(k);
		i = n;
		while ( --i >= 0 ) {
		    j = i;
		    while ( --j >= 0 ) {
			p = nbprm[i][j].getX(); r = nbprm[i][j].getY(); up = 0; down = 0; m = k; 
			while ( --m >= 0 )
			    if ( exists[i][mb=bid[m]] && exists[j][mb] && ( ( sl = ( (im=seqid[i][mb]) > (jm=seqid[j][mb]) ) ? pt2m[mb][im][jm].copy() : pt2m[mb][jm][im].copy() ).getY() > 0 ) ) { pv = Math.min( (pv=nbcdf((s=sl.getX()), p, r*(l=sl.getY()))) , 1-pv+nbpmf(s, p, r*l) ); up += pv*s; down += pv*l; }
			dm[i][j] = up / down;
		    }
		}
		out.write(" " + n); out.newLine(); i = -1; 
		while ( ++i < n ) {
		    out.write((tax.get(i) + BLANK).substring(0,mtl));
		    j = -1; while ( ++j < i ) out.write(" " + df.format(dm[i][j])); out.write(" " + df.format(0)); while ( ++j < n ) out.write(" " + df.format(dm[j][i]));
		    out.newLine(); 
		}
		out.newLine(); 
	    }
	    System.err.println("\r                                                    \r[" + (100*o/bootstrap) + "%] estimating bootstrap pairwise distances ...");
	}
	out.close();
	


    }


    // a class to store pairs of integers
    private static class Point2D {
	private int x, y;
	public Point2D(int x, int y) { this.x = x; this.y = y; }
	public int getX() { return this.x; }
	public int getY() { return this.y; }
	public String toString() { return "(" + this.x + "," + this.y + ")"; }
	public Point2D copy() { return new Point2D(this.x, this.y); }
    }

    // a class to store pairs of doubles
    private static class Point2d {
	private double x, y;
	public Point2d(double x, double y) { this.x = x; this.y = y; }
	public double getX() { return this.x; }
	public double getY() { return this.y; }
	public Point2d copy() { return new Point2d(this.x, this.y); }
    }

    // a class dedicated to the NB parameter estimate
    private static class NBfactory implements Runnable {
	private final static double SMALL = 1E-10;
	private Point2D indexes;
	private Point2d nbp;
	private ArrayList<Point2D> pt2a;
	private double outlierThreshold;
	public NBfactory(Point2D indexes, short[][] seqid, boolean[][] exists, Point2D[][][] pt2m, double outlierThreshold) {
	    this.indexes = indexes.copy(); int i = indexes.getX(), j = indexes.getY(); this.outlierThreshold = outlierThreshold;
	    int k = pt2m.length; pt2a = new ArrayList<Point2D>(k); short im, jm; Point2D sl; int m = k; 
	    while ( --m >= 0 ) if ( exists[i][m] && exists[j][m] && ( ( sl = ((im=seqid[i][m]) > (jm=seqid[j][m])) ? pt2m[m][im][jm] : pt2m[m][jm][im] ).getY() > 0 ) ) pt2a.add( sl.copy() );
	    nbp = new Point2d(0, 0);
	}
	public void run() {
	    boolean ok = false; int s, l, m; double pv, p = 1E-4, r = 1; 
	    while ( ! ok ) {
		ok = true; nbp = MLPDistMT.getNBparameters(pt2a, SMALL, ( r < 0.1 ) ? 10*r : 1 ); 
		if ( (nbp.getX() == 1E-4) && (nbp.getY() == 1) ) { nbp = new Point2d(p, r); break; } // sumk == 0, then return (p,r) from previous iteration
		p = nbp.getX(); r = nbp.getY(); m = pt2a.size(); 
		while ( --m >= 0 ) {
		    //s = pt2a.get(m).getX(); l = pt2a.get(m).getY(); pv = Math.min( (pv=nbcdf(s, p, r*l)) , 1 - pv + nbpmf(s, p, r*l) ); if ( pv < outlierThreshold ) { pt2a.remove(m); ok = false; }
		    if ( Math.min((pv=nbcdf((s=pt2a.get(m).getX()), p, r*(l=pt2a.get(m).getY()))), 1-pv+nbpmf(s, p, r*l)) < outlierThreshold ) { pt2a.remove(m); ok = false; }
		}
	    }
	    pt2a = null;
	}
	public Point2D getIndexes() { return indexes; }
	public Point2d getNBparameters() { return nbp; }
    }



    // computes the square of a double
    static double square(double x) { return x*x; }

    // estimates the p,r parameters of the NB distribution
    static Point2d getNBparameters( ArrayList<Point2D> data ) {
	return getNBparameters(data, SMALL, 1);
    }
 
   // estimates the p,r parameters of the NB distribution
    static Point2d getNBparameters( ArrayList<Point2D> data , double rleft, double rright) {
	double p = 0.5, r = 0.5, lk = Double.NEGATIVE_INFINITY; int size = data.size(); int[] kar = new int[size], lar = new int[size]; double sumk = 0, suml = 0;
	int k, l, m = size; while ( --m >= 0 ) { sumk += (kar[m]=(short)data.get(m).getX()); suml += (lar[m]=(short)data.get(m).getY()); }
	//System.out.println(sumk + " " + suml + " " + (sumk/((double)suml)));
	if ( sumk == 0 ) return new Point2d(1E-4, 1);
	// ##########################################################################################################################################################
	// ######### fast estimate of r with golden section search to minimize lk = sum of -log f(k_m, p, r_m*l_m) with m = 1, ..., k                       #########
	// ######### note that ML maximization leads to p = sumk/(sumk + suml*r)                                                                            #########
	// ######### golden section search starting parameters: r1 ~ 0, r4 ~ 1, r2 = r1 + (r4-r1)/(1+1.6180339887)                                          #########
	// ######### if the golden section search cannot be performed, a naive search is performed, i.e. trying every r = 0 to 1 with step = 0.001          #########
	// ######### golden section search adapted from code.google.com/p/jlabgroovy/source/browse/trunk/GroovyLabSrc/com/nr/min/Golden.java                #########
	// ##########################################################################################################################################################
	/*
	double r1 = rleft, r4 = rright, r2 = r1 + 0.381966*(r4-r1), r3 = r2 + 0.381966*(r2-r1), p1 = sumk/(sumk + suml*r1), p2 = sumk/(sumk + suml*r2), p3 = sumk/(sumk + suml*r3), p4 = sumk/(sumk + suml*r4), lk1 = 0, lk2 = 0, lk3 = 0, lk4 = 0;
	m = size; while ( --m >= 0 ) { lk1 -= lognbpmf((k=kar[m]), p1, r1*(l=lar[m])); lk2 -= lognbpmf(k, p2, r2*l); lk3 -= lognbpmf(k, p3, r3*l); lk4 -= lognbpmf(k, p4, r4*l); }
	*/
	double r1 = rleft, r4 = rright, p1 = sumk/(sumk + suml*r1), p4 = sumk/(sumk + suml*r4), lk1 = 0, lk4 = 0; 
	m = size; while ( --m >= 0 ) { lk1 -= lognbpmf((k=kar[m]), p1, r1*(l=lar[m])); lk4 -= lognbpmf(k, p4, r4*l); }
	double r2 = r4, p2 = p4, lk2 = (lk1 > lk4) ? lk1 : lk4; int c = 20;
	while ( (--c >= 0) && ! ((lk2 < lk1) && (lk2 < lk4)) ) { 
	    if ( (r2 /= 1.6180339) < SMALL ) { r2 *= 1.6180339; break; } p2 = sumk/(sumk + suml*r2); lk2 = 0; m = size; while ( --m >= 0 ) lk2 -= lognbpmf(kar[m], p2, r2*lar[m]);
	}
	double r3 = (lk2 < lk4) ? r2 + 0.381966*(r2-r1) : 1, p3 = sumk/(sumk + suml*r3), lk3 = 0; m = size; while ( --m >= 0 ) lk3 -= lognbpmf(kar[m], p3, r3*lar[m]);
	//System.out.println(r1 + " " + lk1 + "   " + r2 + " " + lk2 + "   " + r4 + " " + lk4);
	if ( (! Double.isInfinite(lk1)) && (! Double.isInfinite(lk2)) && (! Double.isInfinite(lk3)) && (! Double.isInfinite(lk4)) && (! Double.isNaN(lk1)) && (! Double.isNaN(lk2)) && (! Double.isNaN(lk3)) && (! Double.isNaN(lk4)) && (lk2 < lk1) && (lk2 < lk4) ) { // ### golden section search for r
	    while ( Math.abs(r4-r1) > 1E-6 * (Math.abs(r2) + Math.abs(r4)) ) 
		if ( lk3 < lk2 ) { r1 = r2; r2 = r3; r3 = 0.61803399 * r3 + 0.381966 * r4; p2 = p3; p3 = sumk/(sumk + suml*r3); lk2 = lk3; lk3 = 0; m = size; while ( --m >= 0 ) lk3 -= lognbpmf(kar[m], p3, r3*lar[m]); }
		else {             r4 = r3; r3 = r2; r2 = 0.61803399 * r2 + 0.381966 * r1; p3 = p2; p2 = sumk/(sumk + suml*r2); lk3 = lk2; lk2 = 0; m = size; while ( --m >= 0 ) lk2 -= lognbpmf(kar[m], p2, r2*lar[m]); }
	    r = ( lk2 < lk3 ) ? r2 : r3; p = ( lk2 < lk3 ) ? p2 : p3; lk = ( lk2 < lk3 ) ? -lk2 : -lk3; //if ( lk2 < lk3 ) { r = r2; p = p2; lk = -lk2; } else { r = r3; p = p3; lk = -lk3; }
	}
	else { /*System.out.println("naive");*/ // ### naive search for r
	    //r2 = 0; while ( (r2 += 0.001 ) < 1 ) { p2 = sumk / (sumk + suml*r2); lk2 = 0; m = size; while ( --m >= 0 ) lk2 += lognbpmf(kar[m], p2, r2*lar[m]); if ( lk2 > lk ) { lk = lk2; r = r2; p = p2; } } 
	    r = SMALL; p = sumk/(sumk + suml*r); lk = Double.POSITIVE_INFINITY;
	    if ( (! Double.isInfinite(lk1)) && (! Double.isNaN(lk1)) ) { r = r1; p = p1; lk = lk1; }
	    if ( (! Double.isInfinite(lk2)) && (! Double.isNaN(lk2)) && (lk2 < lk) ) { r = r2; p = p2; lk = lk2; }
	    if ( (! Double.isInfinite(lk3)) && (! Double.isNaN(lk3)) && (lk3 < lk) ) { r = r3; p = p3; lk = lk3; }
	    if ( (! Double.isInfinite(lk4)) && (! Double.isNaN(lk4)) && (lk4 < lk) ) { r = r4; p = p4; lk = lk4; }
	    lk = -lk;
	} 
	//System.out.println("# " + size + " " + p + " " + r + " " + lk);
	return new Point2d(p, r);
    }

    // estimates the NB(p,r) PMF, i.e. P(X=x) with X~NB(p,r)
    static double nbpmf(int x, double p, double r) {
	return ( x == 0 ) ? Math.pow(1-p, r) : Math.exp( gammln(r+x) - gammln(x+1) - gammln(r) + x*Math.log(p) + r*Math.log(1-p) );
    }
    
    // estimates the log of the NB(p,r) PMF 
    static double lognbpmf(int x, double p, double r) {
	return ( x == 0 ) ? r * Math.log(1-p) : gammln(r+x) - gammln(x+1) - gammln(r) + x*Math.log(p) + r*Math.log(1-p);
    }

    // estimates the NB(p,r) CDF, i.e. P(X<=x) with X~NB(p,r) 
    static double nbcdf(int x, double p, double r) {
	//int k = x; ++k; double sum = 0; while ( --k >= 0 ) sum += nbpmf(k, p, r); return sum;
	int k = x, k1 = (++k); ++k1; double rk = r + k, cte = r*Math.log(1-p)-gammln(r), logp = Math.log(p), sum = 0; while ( --k >= 0 ) sum += Math.exp( gammln(--rk) - gammln(--k1) + k*logp + cte ); return sum;
    }

    // estimates ln Gamma(x)
    static double gammln(double xx) {
        final double[] LANCZOS = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.001208650973866179, -0.000005395239384953};
        double y_, tmp_, ser_; 
        y_ = xx; tmp_ = xx + 5.5; tmp_ -= (xx + 0.5) * Math.log(tmp_); ser_ = 1.000000000190015;
        int j_ = -1; while ( ++j_ <= 5 ) ser_ += LANCZOS[j_] / ++y_; 
        return -tmp_ + Math.log(2.5066282746310005 * ser_ / xx);
    }
}



