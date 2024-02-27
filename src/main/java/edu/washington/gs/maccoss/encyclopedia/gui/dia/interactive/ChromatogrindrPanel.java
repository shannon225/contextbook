package edu.washington.gs.maccoss.encyclopedia.gui.dia.interactive;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Optional;
import java.util.zip.DataFormatException;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.RowFilter;
import javax.swing.SwingUtilities;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.event.AnnotationChangeListener;
import org.jfree.chart.panel.CrosshairOverlay;
import org.jfree.chart.plot.Crosshair;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.plot.XYPlot;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.TextAnchor;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.parameters.InstrumentSpecificSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.FragmentIonConsistencyCharter;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerProgress;
import edu.washington.gs.maccoss.encyclopedia.gui.massspec.ChromatogramCharter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.CategoricalData;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.EditableXYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BackgroundSubtractionFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import gnu.trove.list.array.TFloatArrayList;

/**
 * click right to approve chromatogram, click left to flag as bad
 * @author searleb
 *
 */
public class ChromatogrindrPanel extends JPanel {
	private static final long serialVersionUID=1L;
	
	private static final float RT_EXTRACTION_MARGIN_IN_SEC=45f;
	
	private final FileChooserPanel diaFileChooser;
	private final FileChooserPanel libraryFileChooser;
	private final JSplitPane mainSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JSplitPane tableSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	
	private final JTable peptideTable;
	private final TableRowSorter<TableModel> rowSorter;
	private final PeptidePrecursorTableModel peptideModel;
	private final JTextField jtfFilter;
	private final JCheckBox jtfNotFilter=new JCheckBox("NOT");
	private final JComboBox<InstrumentSpecificSearchParameters> instrumentCombo=new JComboBox<InstrumentSpecificSearchParameters>(InstrumentSpecificSearchParameters.INSTRUMENTS);

	private final JCheckBox sgSmoothBox;
	private final JCheckBox backgroundSubtractBox;

	private LibraryInterface reference=null;
	private StripeFileInterface dia=null;

	
	public static void main(String[] args) {
		File rawFile=new File("/Users/searleb/Documents/encyclopedia/small_file/bcs_2020jan16_hela_clib_3.mzML");
		File libraryFile=new File("/Users/searleb/Documents/encyclopedia/small_file/pan_human_library.dlib");
		final ChromatogrindrPanel browser=new ChromatogrindrPanel();
		launchBrowserPanel(browser);
		
		browser.updateLibrary(libraryFile);
		browser.updateRaw(rawFile);
		
		browser.pasteTable("KNILLTIGSYK	64.14	2	\n" + "GLGC[+57.021464]SLLFIPLGLVDRR	96.87	3	\n"
				+ "ATTGTQTLLSSGTR	46.81	2	\n" + "STLDIVLANKDR	57.81	2	\n" + "YETLFQALDR	78.86	2	\n"
				+ "AYSFAMGC[+57.021464]WPK	72.86	2	\n" + "YAEEELEQVR	48.26	2	\n"
				+ "RNPDTQWITK	42.77	2	\n" + "RNEFLGELQK	51.95	2	\n"
				+ "ISEEC[+57.021464]IAQWK	57.02	2	\n" + "GGDFSSSDFQSR	44.43	2	\n"
				+ "TVLLLADQMISR	81.11	2	\n" + "LWTSC[+57.021464]APLNIR	69.84	2	\n"
				+ "AGSFITGIDVTSK	66.42	2	\n" + "SMAEDTINAAVK	50.93	2	\n" + "GGGTPDANSLAPPGK	43.07	2	\n"
				+ "ATSLGRPEEEEDELAHR	42.37	3	\n" + "LWDYIDGILIK	89.78	2	\n"
				+ "IPNIYAIGDVVAGPM[+15.994915]LAHK	81.20	3	\n" + "QLFHGTPVTIENFLNWK	88.69	3	\n"
				+ "LPEHC[+57.021464]IEYVR	45.37	2	\n" + "RQESGYLIEEIGDVLLAR	94.88	3	\n"
				+ "DHGLEVLGLVR	72.62	2	\n" + "ISTLTIEEGNLDIQRPK	64.12	3	\n" + "AVQELVHPVVDR	50.54	2	\n"
				+ "MSAYSGITDVIIGMPHR	79.54	3	\n" + "MSQVAPSLSALIGEAVGAR	93.61	3	\n"
				+ "LVTC[+57.021464]TGYHQVR	35.64	2	\n" + "GGPNIITLADIVKDPVSR	84.24	3	\n"
				+ "SNQLFNGHGGHIMPPTQSQFGEMGGK	60.40	4	\n" + "TVEAEAAHGTVTR	31.26	2	\n"
				+ "VAPLWHSSSEVISMADR	68.71	3	\n" + "QVTPDGESDEVGVIPSKR	48.04	3	\n"
				+ "LIETLKPFGVFEEEEELQRR	78.29	4	\n" + "DLEQPSQAAGINLEIIR	80.25	3	\n"
				+ "AYFDLQTHVIQVPQGK	69.42	3	\n" + "AYWLLEEMLTK	94.64	2	\n" + "VPNSVLWLLR	87.95	2	\n"
				+ "AQETGHLVMDVR	46.31	2	\n" + "LAQM[+15.994915]FSDMVLK	67.63	2	\n"
				+ "MLGTEGGEGFVVK	61.23	2	\n" + "QVSGLTVDTEER	46.31	2	\n"
				+ "LVAFGTSHGFGLFDHQQR	69.63	3	\n" + "SGAGEDPPMPASR	38.16	2	\n"
				+ "QVSASELHTSGILGPETLR	64.22	3	\n" + "LVIVDGIAFPFR	92.84	2	\n"
				+ "LIQFC[+57.021464]AIDELGTNYPK	81.94	3	\n" + "SIVEEIEDLVAR	97.14	2	\n"
				+ "HPDSHQLFIGNLPHEVDKSELKDFFQSYGNVVELR	84.82	6	\n"
				+ "DVAHWLGC[+57.021464]SATSTFNFHPNVR	77.60	4	\n" + "LLPALQSTITR	67.46	2	\n"
				+ "DGFVQNVHTPR	44.68	2	\n" + "AEELIQEIQR	59.01	2	\n"
				+ "ELMC[+57.021464]QIEASAK	54.10	2	\n" + "QSVEADINGLR	54.00	2	\n"
				+ "AGLSPANC[+57.021464]QSDR	32.31	2	\n" + "EVWALVQAGIR	78.76	2	\n"
				+ "NFATSLYSMIK	84.59	2	\n" + "IFLYPNAGQLK	69.89	2	\n" + "DLPPVSGSIIWAK	80.74	2	\n"
				+ "MININILSVC[+57.021464]K	80.29	2	\n" + "SIIEC[+57.021464]VDDFR	67.87	2	\n"
				+ "ELPDLEDLMK	84.14	2	\n" + "AIPNQGEILVIR	69.13	2	\n" + "DYIWNTLNSGR	73.05	2	\n"
				+ "GFC[+57.021464]QLVVSSSLR	69.66	2	\n" + "TFQMDDYSLC[+57.021464]GLISHK	77.77	3	\n"
				+ "EYVEPELHINDLWR	78.19	3	\n" + "VKNEGDDFGWGVVVNFSK	79.71	3	\n"
				+ "MLAQPLKDSDVEVYNIIK	75.86	3	\n" + "MVEPQYQELK	50.71	2	\n"
				+ "ATPSENLVPSSAR	46.77	2	\n" + "TESPVLTSSC[+57.021464]R	39.62	2	\n"
				+ "DPNTQSVGNPQR	29.78	2	\n" + "EVMFTEEDVK	54.97	2	\n" + "LTHYDHVLIELTQAGLK	75.81	3	\n"
				+ "EAQELSQNSAIKQDAQSLHGDIPQK	51.97	4	\n" + "EMVSDVDLSFNK	65.81	2	\n"
				+ "VSGVDGYETEGIR	49.14	2	\n" + "ISAFGYLEC[+57.021464]SAK	68.26	2	\n"
				+ "LAPGFDAELIVK	76.12	2	\n" + "NSC[+57.021464]NVGGGGGGFK	32.43	2	\n"
				+ "SRPPEERPPGLPLPPPPPSSSAVFR	69.25	4	\n" + "YDYVLTGYTR	60.16	2	\n"
				+ "DLTPEHLPLLR	71.33	2	\n" + "TSQLLETLNQLSTHTHVVDITR	82.40	4	\n"
				+ "ALLQQQPEDDSK	40.52	2	\n" + "TFSFYLSNIGR	79.80	2	\n" + "ATQELIPIEDFITPLK	96.14	3	\n"
				+ "EAAEAEAEVPVVQYVGER	68.63	3	\n" + "EVMSPLQAMSSYTVAGR	81.15	3	\n"
				+ "IPPLNPGQGPGPNK	49.99	2	\n" + "LVAC[+57.021464]FQGQHGTDAERR	37.40	3	\n"
				+ "GTGAASFDEFGNSK	52.10	2	\n" + "NLNGTLHELLR	66.91	2	\n" + "LESENDEYER	31.44	2	\n"
				+ "AQPTPSSSATQSKPTPVKPNYALK	42.12	4	\n" + "GHDLNEDGLVSWEEYK	65.57	3	\n"
				+ "RLSQIGVENTEENRR	35.68	3	\n" + "APLKPYPVSPSDK	43.30	2	\n" + "NPLVAVYYTNR	62.66	2	\n"
				+ "NLHVVFTMNPSSEGLKDR	61.51	3	\n" + "GAVYSMVEFNGK	63.09	2	\n"
				+ "EYVNSTSEESHDEDEIRPVQQQDLHR	44.33	5	\n" + "EPPADVWTPPAR	58.65	2	\n"
				+ "NADHSMNYQYR	33.96	2	\n" + "NAQEALQAIETK	56.01	2	\n" + "EDLPAENGETK	33.87	2	\n"
				+ "LLEEENQESLR	45.66	2	\n" + "IVFAAGNFWGR	79.69	2	\n" + "LYQGINQLPNVIQALEK	90.46	3	\n"
				+ "QVTSSGVSHGGTVSLQDAVTR	50.02	3	\n" + "THSQGGYGSQGYK	25.90	2	\n"
				+ "MVDENC[+57.021464]VGFDHTVKPVSDMELETPTDKR	60.48	5	\n" + "LAGDPSAGDGAAPR	34.31	2	\n"
				+ "LFAVLEQLSPVR	87.07	2	\n" + "NTLTNIAM[+15.994915]RPGLEGYALPR	68.67	3	\n"
				+ "KDPGVPNSAPFK	43.93	2	\n" + "SEIC[+57.021464]TEEPQK	32.30	2	\n"
				+ "C[+57.021464]DPAGYYC[+57.021464]GFK	53.99	2	\n" + "ALEQQVEEMK	47.84	2	\n"
				+ "NLTGDVC[+57.021464]AVMR	60.78	2	\n" + "TC[+57.021464]LIC[+57.021464]ADTFR	62.03	2	\n"
				+ "LAQDGAHVVVSSR	36.45	2	\n" + "LVSWYTLMEGQEPIAR	87.46	3	\n" + "EDLYLKPIQR	52.89	2	\n"
				+ "VQFAPEKPGPQPSAETTR	44.37	3	\n" + "EQSGTIYLQHADEEREK	39.95	3	\n"
				+ "TDAEATDTEATET	32.20	2	\n" + "QDLPALEEKPR	46.14	2	\n" + "LVPGGGATEIELAK	59.74	2	\n"
				+ "VLHMVGDKPVFSFQPR	64.95	3	\n" + "WEEVQSYIR	59.78	2	\n" + "VFEVNASNLEK	56.81	2	\n"
				+ "LTGTIQNDILK	59.81	2	\n" + "SVSGTDVQEEC[+57.021464]R	31.45	2	\n"
				+ "LPPNTNDEVDEDPTGNK	42.33	3	\n" + "NVALSGVLEVVR	79.42	2	\n" + "QLQLAQEAAQK	44.98	2	\n"
				+ "TIEDLDENQLKDEFFK	74.48	3	\n" + "EVWDYVFFK	87.08	2	\n"
				+ "ISM[+15.994915]PDLDLNLK	67.16	2	\n" + "WDLSAQQIEER	62.26	2	\n"
				+ "TASSVIELTC[+57.021464]TK	55.26	2	\n" + "QRPGQQVATC[+57.021464]VR	31.39	2	\n"
				+ "VFSDEVQQQAQLSTIR	61.03	3	\n" + "RLEFPSGETIVMHNPK	58.85	3	\n"
				+ "FNC[+57.021464]EENQHSDSC[+57.021464]YK	29.61	3	\n" + "VNPYEEVDQEK	45.28	2	\n"
				+ "SLEDALAEAQR	58.97	2	\n" + "C[+57.021464]HWSDMFTGR	58.77	2	\n"
				+ "QLDTVNFFLK	80.84	2	\n" + "QQPDTEAVLNGK	44.00	2	\n"
				+ "ATSITVTGSGSC[+57.021464]R	37.35	2	\n" + "VSNQVAVNMYK	47.48	2	\n"
				+ "QLGELLTDGVR	65.94	2	\n" + "EADASPASAGIC[+57.021464]R	37.79	2	\n"
				+ "TPLHMAASEGHASIVEVLLK	74.64	3	\n");
		
		/*browser.pasteTable("ELNVMFIETSAK	75.89	2	\n" + "QHYVLAGASGSPGEEVAIRPSTAPR	53.64	4	\n"
				+ "ALLPILQWHK	81.15	2	\n" + "IATEFNQLQFHAVQSK	62.84	3	\n"
				+ "EQM[+15.994915]QPTHPIR	24.31	2	\n" + "LDNTTAAVQELGR	53.34	2	\n"
				+ "EDGTIFHPHSGLC[+57.021464]LSAYR	55.27	3	\n" + "TVTPAMVEGIYK	73.42	2	\n"
				+ "IPGIYVLSLEIGK	89.67	2	\n" + "LTAYAM[+15.994915]TIPFVR	74.86	2	\n"
				+ "LTSPVINTSLDTK	59.47	2	\n" + "LDLTGTSGTAVPAR	53.99	2	\n" + "DYTEGWVEFR	73.85	2	\n"
				+ "IAALENELTFLR	83.55	2	\n" + "DFQYNEEEMK	49.04	2	\n" + "SPPSAGYLVMVSR	67.14	2	\n"
				+ "LEELELDEQQK	53.37	2	\n" + "APGFAHLAGLDK	57.70	2	\n"
				+ "KVVLMQC[+57.021464]NIESVEEGVK	62.32	3	\n" + "LSKDQFALAMYFIQQK	88.02	3	\n"
				+ "MATYLTGELTATSEDYK	69.09	3	\n" + "QLC[+57.021464]EMEAC[+57.021464]R	38.90	2	\n"
				+ "RIIDDSEITKEDDALWPPPDR	66.88	4	\n" + "LSSLGGVNSLGVSSLEHITHSLLGR	90.05	4	\n"
				+ "HVAMTLLDTEQSYVESLR	78.91	3	\n" + "TSSEPEFNSLPR	54.43	2	\n" + "LQDEGQEAEGEK	26.93	2	\n"
				+ "NVTEMAMNPHIK	49.05	2	\n" + "EQHGLQLQSEINQLHSK	54.39	3	\n"
				+ "QQLPTFLQQMQNPDTLSAMSNPR	74.63	4	\n" + "SIIQQHNLETLENDIK	66.66	3	\n"
				+ "HEM[+15.994915]LPEFYK	51.64	2	\n" + "KTDFSWEEERNFGASLLLPGLK	84.37	4	\n"
				+ "PAMLHLPSEQGAPETLQR	61.66	3	\n" + "GGTFQM[+15.994915]GGSSSHNRPSGSNVDTLLR	47.76	4	\n"
				+ "EQGLRDIASTPHELYR	60.82	3	\n" + "SC[+57.021464]C[+57.021464]TFC[+57.021464]GVLR	56.53	2	\n"
				+ "GQVTGALLLSVVGGK	83.71	2	\n" + "MHIAQDINQDNLQLFLNSYNGRR	76.48	4	\n"
				+ "IHQQELEVGISSHQPSFAALNR	61.81	4	\n" + "KNENLQNLLC[+57.021464]GSGAGVISK	70.50	3	\n"
				+ "TYVSPTHVGSGAYGSVC[+57.021464]SAIDKR	58.51	4	\n" + "VGDAIPAVEVFEGEPGNKVNLAELFK	92.47	4	\n"
				+ "FNAC[+57.021464]FESVATNIDEIYK	82.51	3	\n" + "LASTNSSVLGADLPSSMKEK	59.06	3	\n"
				+ "HKPGIVQETTFDLGGDIHSGTALPTSK	70.41	4	\n" + "FGTFPGNYVAPV	83.12	2	\n"
				+ "SSYIVSQIAVAYHNIR	75.14	3	\n" + "YQLEIKIPETYPFNPPK	83.52	3	\n"
				+ "GRANHSAFLFGFGDGGGGPTQTMLDR	80.74	4	\n" + "TNSQLDTSIQR	39.53	2	\n"
				+ "TGFTPLHIAAHYENLNVAQLLLNR	88.35	4	\n"
				+ "QPTFC[+57.021464]SHC[+57.021464]TDFIWGFGK	86.85	3	\n"
				+ "RPFAVTTQSFGSNAEGQHSGFGPQPNPEK	56.71	5	\n"
				+ "AAVHYTVGC[+57.021464]LC[+57.021464]EEVALDK	63.54	3	\n" + "DIKPENLLISHNDVLK	72.72	3	\n"
				+ "DSQMQNPYSR	36.72	2	\n" + "RPQAVIEDAVATSGVSTLSSTVSHDSQSAHR	64.52	5	\n"
				+ "LVRPEVDVMC[+57.021464]TAFHDNEETFLK	78.27	4	\n" + "QNLEPLFDSYTSELRR	70.25	3	\n"
				+ "VSSVPNTSQSYAK	36.36	2	\n" + "AEYINFLENLK	83.85	2	\n"
				+ "QMSVKEDLDKVEPAVIEAQNAVK	73.54	4	\n" + "LDAYKADDPTMGEGPDK	44.20	3	\n"
				+ "IHNEMASTSDK	26.19	2	\n" + "TVLMLADQM[+15.994915]ISR	66.40	2	\n"
				+ "GQVLPAHTLLNTVDVELIYEGVK	93.21	4	\n" + "KFEEGSFANSTDQEPTRPQPGGGDVR	49.02	4	\n"
				+ "SSESEFTQYTTHHILK	53.73	3	\n" + "EQLSLPAEFPDK	69.29	2	\n"
				+ "GGLGGGYGGASGM[+15.994915]GGITAVTVNQSLLSPLVLEVDPNIQAVR	92.40	6	\n"
				+ "EWYIGYYQGR	68.72	2	\n" + "FNAESQGC[+57.021464]NHEEDAGVR	33.28	3	\n"
				+ "SANFLDHLYVGIPRPSGEK	74.71	3	\n" + "SVYHQLFMSSLLMDLK	92.46	3	\n"
				+ "RNPAGSVVMER	36.03	2	\n" + "LVLVSPTSEQYDSLLR	80.28	3	\n"
				+ "SLC[+57.021464]PETWPTWAGRPQDGVAVLVR	83.32	4	\n"
				+ "TAASGIPYHSEVPVSLKEAVC[+57.021464]EVALDYKK	81.91	5	\n"
				+ "GTDLWLGVDALGLNIYEKDDKLTPK	87.57	4	\n" + "IVLPGNFLYC[+57.021464]TFYGR	93.70	3	\n"
				+ "DDKESVPISDTIIPAVPPPTDLR	75.96	4	\n" + "LGTVYC[+57.021464]QASFPGANIIGNK	74.88	3	\n"
				+ "SEEEQSSSSVKKDETNVK	25.71	3	\n" + "TGTYRQLFHPEQLITGKEDAANNYAR	61.93	5	\n"
				+ "KTGVAGEDMQDNSGTYGK	33.97	3	\n" + "KVIDQQNGLYR	38.14	2	\n" + "SASLVVPSDIPK	53.32	2	\n"
				+ "QC[+57.021464]C[+57.021464]VLFDFVSDPLSDLK	97.64	3	\n"
				+ "LAQLEEAKQASIQHIQNAIDTEK	68.53	4	\n" + "MVFFVQNEPPHQIFK	74.71	3	\n"
				+ "LEVERDNLAQDLATVR	73.22	3	\n" + "VALGNTWKENLTELSGGQR	68.33	3	\n"
				+ "DVFGTNQLVGC[+57.021464]R	63.35	2	\n" + "DVQM[+15.994915]LQDAISK	55.22	2	\n"
				+ "ENALQDSILAR	61.71	2	\n" + "DHPFGFVAVPTKNPDGTMNLMNWEC[+57.021464]AIPGKK	82.83	5	\n"
				+ "LEQC[+57.021464]PLQLNNPFNEYSK	76.26	3	\n"
				+ "LSQERPGVLLNQFPC[+57.021464]ENLLTVK	78.55	4	\n" + "KGLPDQELFSLNEGVR	75.20	3	\n"
				+ "QDPGDNWEEGGGGGGGMEK	42.15	3	\n" + "EALELTDTGLLSGSEER	73.99	3	\n"
				+ "NQIKVDLVDENFTELR	73.35	3	\n" + "YHDSDEATAAR	23.83	2	\n" + "EWEEAELQAK	50.77	2	\n"
				+ "VC[+57.021464]LYPGFVDVK	70.15	2	\n" + "FRPLQLETINVTM[+15.994915]AGK	68.68	3	\n"
				+ "AELQAQLAALSTK	66.30	2	\n" + "EAFLVNSDLTLR	75.89	2	\n" + "EALQSDLLEMK	67.72	2	\n"
				+ "C[+57.021464]TC[+57.021464]GFSAIMNR	58.74	2	\n" + "VMPIC[+57.021464]LPSKDYAEVGR	59.91	3	\n"
				+ "STPYSAYDPETYTGHWK	61.76	3	\n" + "ETVVISPPC[+57.021464]TGSSEHWKPELEEK	60.30	4	\n"
				+ "VRFLEQQNAALAAEVNR	62.25	3	\n" + "LC[+57.021464]DSGELVAIKK	47.99	2	\n"
				+ "AFENLLGQALTK	78.65	2	\n" + "VC[+57.021464]DAC[+57.021464]FNDLQG	57.96	2	\n"
				+ "NLAATLQDIETK	68.54	2	\n" + "KC[+57.021464]NLVPTDEITVYYK	64.79	3	\n"
				+ "SPNNFLSYYR	61.50	2	\n" + "DSVASTITGVMDK	65.86	2	\n" + "ELDREAQAEYLLQVR	63.42	3	\n"
				+ "GVPNVISEDTLK	61.50	2	\n" + "KLAPEEC[+57.021464]FSPLDLFNK	83.28	3	\n"
				+ "LVIVSLMELFK	104.72	2	\n" + "MQYAPNTQVEILPQGHESPIFK	72.76	4	\n"
				+ "NSAEAIIHGLSSLTAC[+57.021464]QLR	82.76	3	\n" + "RLSEDYGVLKTDEGIAYR	57.45	3	\n"
				+ "KIPVFHNGSTPTLGETPK	51.84	3	\n" + "EPVC[+57.021464]AALNSAILESQNLPK	80.59	3	\n"
				+ "NIFHLFHDVVPTYHK	72.07	3	\n" + "KVADALTNAVAHVDDMPNALSALSDLHAHK	83.90	5	\n"
				+ "TTANLAVDVIASSFGQTR	86.77	3	\n" + "DENSQLVAIVLR	82.51	2	\n"
				+ "DNVGEEVDAEQLIQEAC[+57.021464]R	76.20	3	\n" + "HIAEDADRKYEEVAR	35.87	3	\n"
				+ "LVHTNEVTVLLGDNWFAK	84.85	3	\n" + "WC[+57.021464]FLDATTASR	72.39	2	\n"
				+ "SFVHPKPGAAGSVGAGLIPISSELC[+57.021464]YR	73.89	4	\n" + "GKLDGNQDLIR	39.72	2	\n"
				+ "DGTGVVEFVRK	51.73	2	\n" + "KVQGGALEDSQLVAGVAFKK	60.86	3	\n"
				+ "TVIEQQPVLC[+57.021464]EVFC[+57.021464]R	79.69	3	\n" + "NRLLPQGLAVYASPENK	61.12	3	\n"
				+ "AC[+57.021464]PHMATC[+57.021464]GNVLFEGR	55.49	3	\n" + "MVSGM[+15.994915]YLGELVR	79.25	2	\n"
				+ "AFYNNVLGEYEEYITK	82.26	3	\n" + "HSMLFIEASAK	57.02	2	\n"
				+ "VLPVYMNC[+57.021464]LLK	81.17	2	\n" + "DVGAQILLHSHK	50.35	2	\n"
				+ "GQTC[+57.021464]VVHYTGMLEDGK	54.03	3	\n" + "THGTC[+57.021464]AENFYR	36.70	2	\n"
				+ "KAHGLLAEENR	31.28	2	\n" + "SNIVTSINFSK	61.40	2	\n"
				+ "DLKPENLLC[+57.021464]M[+15.994915]GPELVK	72.48	3	\n"
				+ "C[+57.021464]AQSAYC[+57.021464]NTK	25.65	2	\n" + "GSSYLGIPFNPSK	72.52	2	\n"
				+ "RVDFHDVQDYADNIK	59.26	3	\n" + "EGEETLRIEDILEVIEK	93.97	3	\n"
				+ "ILC[+57.021464]HMQLSSAQVEQLR	60.10	3	\n" + "IM[+15.994915]GLDLPDGGHLTHGYMSDVKR	60.87	4	\n"
				+ "NMINTFVPSGK	63.45	2	\n" + "SISLLC[+57.021464]LEGLQK	81.11	2	\n"
				+ "NTFYETLPVAINGNGPTK	73.16	3	\n" + "EAEAAFLNVYK	65.82	2	\n" + "TATATLMLQNR	52.11	2	\n"
				+ "SVAC[+57.021464]DVGYPALK	55.60	2	\n" + "DLFNVDAFKLESLEAK	85.86	3	\n"
				+ "LQAFGNEC[+57.021464]SIEQMEHVR	62.48	3	\n" + "VC[+57.021464]IEHHTFFR	46.96	2	\n"
				+ "LQEMEILYKK	55.75	2	\n" + "LQGC[+57.021464]VSVQVNAGPLAYAR	64.61	3	\n"
				+ "AVWDAFC[+57.021464]ANR	67.62	2	\n" + "SNSLSEQLAINTSPDAVK	60.13	3	\n"
				+ "KAYWQVHLDQVEVASGLTLC[+57.021464]K	75.61	4	\n" + "YSISLSPPEQQK	57.70	2	\n"
				+ "HLQPSQAQPETSIFDVLK	75.82	3	\n" + "RAQAATWANDGLDAEPSK	48.16	3	\n"
				+ "LLEVTADLAER	65.60	2	\n" + "KGPSFADMEVLYWTHVK	80.05	3	\n" + "HTAFATFPNEK	47.26	2	\n"
				+ "TTHQDEEVFK	33.68	2	\n" + "VTLILELLQHK	85.50	2	\n" + "GPSFADMEVLYWTHVK	85.61	3	\n"
				+ "EASADLSPYVR	52.23	2	\n" + "EQLELFQNIRPLFINK	86.97	3	\n" + "LTDTTFLPSSK	58.35	2	\n"
				+ "SIHQIRPSC[+57.021464]AFPVC[+57.021464]HDTEER	46.07	4	\n"
				+ "EQC[+57.021464]DFSNSLK	43.67	2	\n" + "TAWVFDDKYKRPGYGAYDAFK	66.81	4	\n"
				+ "IPM[+15.994915]PVNFNEPLSMLQR	79.76	3	\n" + "IC[+57.021464]LSISGHHPETWQPSWSIR	78.29	4	\n"
				+ "GIEDDLMDLIK	91.68	2	\n" + "ELIQKELTIGSK	53.28	2	\n"
				+ "IQDLKPQC[+57.021464]VVFLNIPR	78.15	3	\n" + "TTVTQSVADSLK	53.67	2	\n"
				+ "LQLQGLDLSSR	68.34	2	\n" + "MEEGGNLGGLIK	58.57	2	\n"
				+ "SC[+57.021464]LLHQFIEK	58.05	2	\n" + "LPPYSAGDGAELSTPGGKLPR	58.63	3	\n"
				+ "TPDSFEESQGEEIGKVER	49.08	3	\n" + "SGPIFIVVPNGK	73.97	2	\n"
				+ "YYYDGDMIC[+57.021464]K	53.86	2	\n" + "LADGGATNQGRVEIFYR	59.76	3	\n"
				+ "LREMLIC[+57.021464]TNMEDLREQTHTR	64.61	4	\n" + "KQLAEQEELER	35.71	2	\n"
				+ "DIPPILRPSLHSETWEIPFEK	75.87	4	\n" + "LQFLAGC[+57.021464]FGLGTVGHTGGK	74.35	3	\n"
				+ "YGLLPSHASYL	72.56	2	\n" + "MRGEAEAFAIGAR	50.76	2	\n" + "YWEMMPPTILIDLLKK	104.06	3	\n"
				+ "TTLADC[+57.021464]LISSNGIISSR	77.79	3	\n" + "VLSGLGGAAASSHR	40.19	2	\n"
				+ "SGMYTVAMAYC[+57.021464]GSGNNK	59.50	3	\n" + "RDLMAC[+57.021464]AQTGSGK	34.77	2	\n"
				+ "ATFSPIVTVEPR	68.07	2	\n" + "M[+15.994915]GITEYNNQC[+57.021464]R	40.36	2	\n"
				+ "YTSQLPPLTAFILPSGGK	93.52	3	\n" + "IC[+57.021464]TGQVPSAEDEPAPKK	39.62	3	\n"
				+ "IYEGAYHVLHK	42.78	2	\n" + "TVIVNM[+15.994915]VDVAK	51.89	2	\n"
				+ "DTPGC[+57.021464]ATTPPHSQASSVR	33.18	3	\n" + "ITHYNYLILSK	63.38	2	\n"
				+ "EGTEAEPLPLR	54.72	2	\n" + "EKPPGASVELVEYLESR	82.25	3	\n"
				+ "DIQNTQC[+57.021464]LLNVEHLSAGC[+57.021464]PHVTLQFADSK	80.70	5	\n"
				+ "LRSEMIEAIR	51.05	2	\n" + "DAVEKPQEFTIVAFVK	79.07	3	\n" + "DWQSYYYHHPQDRDR	44.97	3	\n"
				+ "MHEDINEEWISDKTR	51.27	3	\n" + "KLTAGEAC[+57.021464]AQGLVTEVFPDSTFQK	83.32	4	\n"
				+ "HLIPAANTGESKVFYYK	57.81	3	\n" + "YRWVEQHLGPQFVER	66.26	3	\n"
				+ "LVHPGVAEVVFVK	64.35	2	\n" + "EVDALDGLC[+57.021464]SR	57.14	2	\n"
				+ "LTRDDVIQIC[+57.021464]GPADGIR	61.97	3	\n" + "INLAAATHSAPPFPAAVGSQR	70.88	3	\n"
				+ "RLDDSLLYLR	68.06	2	\n" + "VILEDVAMLHIKPDQFTYTSDHFETIMK	83.98	5	\n"
				+ "QAIKELPQFATGENLPR	64.80	3	\n" + "NAGPIANYLQQVMQEAR	90.15	3	\n"
				+ "LFFFMAPPHQLEFIQK	89.84	3	\n" + "AIGVGLGFELQR	74.40	2	\n"
				+ "VSGLM[+15.994915]MANHTSISSLFER	76.91	3	\n"
				+ "DVLHQNFESYKPEVQELIC[+57.021464]VADR	82.89	4	\n"
				+ "QLILVGDHC[+57.021464]QLGPVVM[+15.994915]C[+57.021464]K	66.93	3	\n"
				+ "FHFFEDQLR	65.01	2	\n" + "MILEIQSMQGK	57.00	2	\n" + "ILMAIDSELVDR	76.17	2	\n"
				+ "THNVHVEIEQR	32.61	2	\n" + "QISRPSAAGINLM[+15.994915]IGSTR	59.60	3	\n"
				+ "EGTFQGLISLR	76.30	2	\n" + "LAAQPLC[+57.021464]MTQPTASGTLR	63.09	3	\n"
				+ "ERFQFPAQVTDVSENAK	64.65	3	\n" + "HPALSPVYLGLLTDWGQR	89.66	3	\n"
				+ "FNHEQHEYYHTHIPNIFQK	55.81	4	\n" + "TLVEQLLSLLNSSPGPPTR	105.20	3	\n"
				+ "NLLQLC[+57.021464]PQSLEALAVR	89.89	3	\n" + "FGVEQDVDMVFASFIRK	95.97	3	\n"
				+ "SHC[+57.021464]IAEVENDEM[+15.994915]PADLPSLAADFVESKDVC[+57.021464]K	73.80	5	\n"
				+ "ITVVGVGQVGM[+15.994915]AC[+57.021464]AISILGK	84.22	3	\n" + "VWDLQAALDPR	79.26	2	\n"
				+ "GTELWLGVDALGLNIYEHDDKLTPK	86.26	4	\n" + "ISIVENC[+57.021464]FGAAGQPLTIPGR	79.13	3	\n"
				+ "VPIWDQDIQFLPGSQK	85.15	3	\n" + "VESTDVSDLLHQYREANQ	67.86	3	\n"
				+ "NHEEEM[+15.994915]KDLR	23.09	2	\n" + "PSLSHLLSQYYGAGVAR	78.44	3	\n"
				+ "LTPVSLSNSPIK	62.43	2	\n" + "ILEDVVGVPEK	61.62	2	\n"
				+ "GDFFPPERPQQLPHGLGGIGMGLGPGGQPIDANHLNK	76.21	6	\n" + "VEGTDGHEAFLLTEGSEEK	57.66	3	\n"
				+ "KINESTQNWHQLENIGNFIK	76.36	4	\n" + "VLDLIVNGISINSAYTSK	85.03	3	\n"
				+ "WLISTDLDQPAAIAVNPK	83.56	3	\n"
				);
		*/
	}

	public static void launchBrowserPanel(final ChromatogrindrPanel browser) {
		final JFrame dialog=new JFrame("Chromatogrindr Browser");

		JMenuBar bar=new JMenuBar();
		JMenu fileMenu=new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		bar.add(fileMenu);
		
		JMenuItem openElib=new JMenuItem("Open Raw Data File...");
		openElib.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				browser.askForData();
			}
		});
		openElib.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(openElib);
		
		dialog.setJMenuBar(bar);
		
		dialog.getContentPane().add(browser, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(1900, 1030);
		dialog.setVisible(true);

		dialog.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				browser.copyTable();
				JOptionPane.showMessageDialog(dialog, "Data copied, remember to paste into spreadsheet!");
			}
		});
	}

	public ChromatogrindrPanel() {
		super(new BorderLayout());
		
		peptideModel=new PeptidePrecursorTableModel();
		peptideTable=new JTable(peptideModel) {
			private static final long serialVersionUID = 1L;

			@Override
			public Object getValueAt(int row, int column) {
				if (column==0) return row+1;
				return super.getValueAt(row, column);
			}
		};
		rowSorter=new TableRowSorter<TableModel>(peptideTable.getModel());
		peptideTable.setRowSorter(rowSorter);

		jtfFilter=new JTextField();
		jtfFilter.getDocument().addDocumentListener(new DocumentListener() {
			@Override
			public void insertUpdate(DocumentEvent e) {
				updateFilter();
			}

			@Override
			public void removeUpdate(DocumentEvent e) {
				updateFilter();
			}

			@Override
			public void changedUpdate(DocumentEvent e) {
				throw new UnsupportedOperationException("Not supported yet.");
			}
		});
		jtfNotFilter.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				updateFilter();
			}
		});

		peptideTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent e) {
				updateToSelectedPeptide();
			}
		});

		libraryFileChooser=new FileChooserPanel(null, "Reference", new SimpleFilenameFilter(LibraryFile.DLIB, LibraryFile.ELIB), false) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filenames) {
				super.update(filenames);
				if (filenames!=null&&filenames.length>0&&filenames[0]!=null) {
					updateLibrary(filenames[0]);
				}
			}
		};

		diaFileChooser=new FileChooserPanel(null, "Raw file", StripeFileGenerator.getFilenameFilter(), true) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filenames) {
				super.update(filenames);
				if (filenames!=null&&filenames.length>0&&filenames[0]!=null) {
					updateRaw(filenames[0]);
				}
			}
		};
		

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(diaFileChooser);
		options.add(libraryFileChooser);
		options.add(instrumentCombo);
		
		instrumentCombo.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				updateToSelectedPeptide();
			}
		});

		sgSmoothBox=new JCheckBox("Smooth");
		sgSmoothBox.setSelected(true);
		sgSmoothBox.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				updateToSelectedPeptide();
			}
		});
		backgroundSubtractBox=new JCheckBox("Background Subtract");
		backgroundSubtractBox.setSelected(true);
		backgroundSubtractBox.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				updateToSelectedPeptide();
			}
		});
		options.add(sgSmoothBox);
		options.add(backgroundSubtractBox);
		
		JPanel buttons=new JPanel(new FlowLayout());
		options.add(buttons);
		
		JButton copyButton = new JButton("Copy");
		copyButton.addActionListener(new ActionListener() {	
			@Override
			public void actionPerformed(ActionEvent e) {
				copyTable();
			}
		});
		buttons.add(copyButton);
		
		JButton pasteButton=new JButton("Paste");
		pasteButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
				try {
					String clip=(String)clipboard.getData(DataFlavor.stringFlavor);
					pasteTable(clip);
					
				} catch (IOException | UnsupportedFlavorException ex) {
					Logger.errorLine("Error reading clipboard!");
					Logger.errorException(ex);
				}
			}
		});
		buttons.add(pasteButton);
		
		peptideTable.addKeyListener(new KeyListener() {
			
			@Override
			public void keyTyped(KeyEvent e) {
			}
			
			@Override
			public void keyReleased(KeyEvent e) {
			}
			
			@Override
			public void keyPressed(KeyEvent e) {
				// click right to approve chromatogram, click left to flag as bad

				boolean hasAction=false;
				if (e.getKeyCode()==KeyEvent.VK_LEFT) {
					hasAction=true;
				} else if (e.getKeyCode()==KeyEvent.VK_RIGHT) {
					hasAction=true;
				} else if (e.getKeyCode()==KeyEvent.VK_ESCAPE) {
					hasAction=true;
				}
				
				if (!hasAction) return;

				int rowIndex = peptideTable.getRowSorter().convertRowIndexToModel(peptideTable.getSelectedRow());
				InteractivePeptidePrecursor peptide=peptideModel.getSelectedRow(rowIndex);
				
				boolean increment=true;
				if (e.getKeyCode()==KeyEvent.VK_LEFT) {
					peptide.setIsPassing(false);
				} else if (e.getKeyCode()==KeyEvent.VK_RIGHT) {
					peptide.setIsPassing(true);
				} else if (e.getKeyCode()==KeyEvent.VK_ESCAPE) {
					peptide.removeIsPassing();
					peptide.setRtRangeInSecs(null);
					updateToSelectedPeptide();
					increment=false;
				}
				peptideModel.fireTableRowsUpdated(rowIndex, rowIndex);
				
				if (increment) {
					int nextRow=peptideTable.getSelectedRow()+1;
					if (nextRow<peptideTable.getRowCount()) {
						peptideTable.setRowSelectionInterval(nextRow, nextRow);
					}
				}
			}
		});
		
		JPanel searchPanel=new JPanel(new BorderLayout());
		searchPanel.add(new JLabel("Search:"), BorderLayout.WEST);
		searchPanel.add(jtfFilter, BorderLayout.CENTER);
		searchPanel.add(jtfNotFilter, BorderLayout.EAST);

		JPanel left=new JPanel(new BorderLayout());
		left.add(options, BorderLayout.NORTH);
		left.add(new JScrollPane(peptideTable, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED), BorderLayout.CENTER);
		left.add(searchPanel, BorderLayout.SOUTH);
		
		mainSplit.setLeftComponent(left);
		mainSplit.setRightComponent(new JLabel("Select a peptide!"));
		
		setLayout(new BorderLayout());
		add(mainSplit, BorderLayout.CENTER);
	}

	public void askForData() {
		diaFileChooser.askForFiles();
	}

	private void updateFilter() {
		String text=jtfFilter.getText();

		if (text.trim().length()==0) {
			rowSorter.setRowFilter(null);
		} else if (jtfNotFilter.isSelected()) {
			rowSorter.setRowFilter(RowFilter.notFilter(RowFilter.regexFilter("(?i)"+text)));
		} else {
			rowSorter.setRowFilter(RowFilter.regexFilter("(?i)"+text));
		}
	}

	public void pasteTable(String clip) {
		peptideModel.paste(clip);
		if (peptideTable.getRowCount()>0) {
			peptideTable.setRowSelectionInterval(0, 0);
		}
		peptideTable.requestFocus();
	}

	private void copyTable() {
		String copyString=peptideModel.copy();
		StringSelection stringSelection = new StringSelection(copyString);
		Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
		clipboard.setContents(stringSelection, null);
	}
	
	public void updateLibrary(final File f) {
		SwingWorkerProgress<LibraryFile> worker=new SwingWorkerProgress<LibraryFile>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Library") {
			@Override
			protected LibraryFile doInBackgroundForReal() throws Exception {
				LibraryFile.OPEN_IN_PLACE=true;
				LibraryInterface ilib=BlibToLibraryConverter.getFile(f);
				LibraryFile.OPEN_IN_PLACE=false;
				if (!(ilib instanceof LibraryFile)) {
					throw new EncyclopediaException("Sorry, can't load this type of library file "+ilib.getClass().getName());
				}
				LibraryFile library=(LibraryFile)ilib;
				return library;
			}
			@Override
			protected void doneForReal(LibraryFile t) {
				Logger.logLine("Finished loading library, updating GUI");
				reference=t;
				updateToSelectedPeptide();
			}
		};
		worker.execute();
	}

	public void updateRaw(final File f) {
		SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Raw File") {
			@Override
			protected Nothing doInBackgroundForReal() throws Exception {
				
				dia=StripeFileGenerator.getFile(f, getParameters(), true);

				Logger.logLine("Read "+dia.getOriginalFileName()+", ("+dia.getRanges().size()+" total windows)");
				return Nothing.NOTHING;
			}
			@Override
			protected void doneForReal(Nothing t) {
				updateToSelectedPeptide();
			}
		};
		worker.execute();
	}

	private SearchParameters getParameters() {
		return instrumentCombo.getItemAt(instrumentCombo.getSelectedIndex()).getDefaultParameters();
	}
	
	public void updateToSelectedPeptide() {
		int[] selection=peptideTable.getSelectedRows();
		if (selection.length<=0) return;
		
		InteractivePeptidePrecursor entry=peptideModel.getSelectedRow(peptideTable.convertRowIndexToModel(selection[0]));
		resetPeptide(entry);
	}
	public void resetPeptide(final InteractivePeptidePrecursor entry) {
		int location=mainSplit.getDividerLocation();
		if (location<=5) {
			location=200;
		}
		SearchParameters parameters=getParameters();
		
		JPanel dataPanel=new JPanel(new GridLayout(0, 2));
		dataPanel.setBackground(Color.WHITE);
		
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] primaryIonObjects=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false);
		
		loadDataBlock: try {
			
			ArrayList<XYTrace> fragmentTraces=new ArrayList<>();
			ArrayList<XYTrace> precursorTraces=new ArrayList<>();
			
			float rtInSec=entry.getRetentionTimeInSec();
			float minRTInSec = rtInSec-RT_EXTRACTION_MARGIN_IN_SEC;
			float maxRTInSec = rtInSec+RT_EXTRACTION_MARGIN_IN_SEC;

			Range rtRange=entry.getRTRange();
			if (rtRange!=null&&rtRange.getRange()>0.0f) {
				minRTInSec=rtRange.getStart();
				maxRTInSec=rtRange.getStop();
			}
			Logger.logLine("Graphing "+entry.getPeptideModSeq()+" ("+primaryIonObjects.length+"), ["+minRTInSec+" to "+maxRTInSec+"]...");
			
			// get precursor traces
			ArrayList<PrecursorScan> precursors=dia.getPrecursors(minRTInSec, maxRTInSec);
			Collections.sort(precursors);
			ArrayList<PrecursorScan> trimmedPrecursors=new ArrayList<>();
			for (PrecursorScan spectrum : precursors) {
				if (entry.getPrecursorMZ()>spectrum.getIsolationWindowLower()&&entry.getPrecursorMZ()<spectrum.getIsolationWindowUpper()) {
					trimmedPrecursors.add(spectrum);
				}
			}
			precursors=trimmedPrecursors;
			XYTraceInterface[] traceArray=ChromatogramExtractor.extractPrecursorChromatograms(parameters.getPrecursorTolerance(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), precursors, true, false);
			for (int i = 0; i < traceArray.length; i++) {
				if (traceArray[i] instanceof XYTrace) {
					precursorTraces.add((XYTrace)traceArray[i]);
				}
			}
			
			// get fragment traces
			ArrayList<FragmentScan> scans=dia.getStripes(entry.getPrecursorMZ(), minRTInSec, maxRTInSec, false);
			if (scans.size()==0) {
				dataPanel.add(new JLabel("no MSMS found from ["+minRTInSec+" to "+maxRTInSec+"] for "+entry.getPrecursorMZ()+" m/z!"));
				break loadDataBlock;
			}
			Collections.sort(scans);
			double[][] allMasses=new double[scans.size()][];
			float[][] allDeltaMasses=new float[scans.size()][];
			float[][] allIntensities=new float[scans.size()][];
			float[] retentionTimes=new float[scans.size()];
			for (int i=0; i<scans.size(); i++) {
				FragmentScan scan=scans.get(i);
				Triplet<double[], float[], float[]> results=extract(scan, primaryIonObjects, parameters);
				double[] masses=results.x;
				float[] deltaMasses=results.y;
				float[] intensities=results.z;
				
				allMasses[i]=masses;
				allDeltaMasses[i]=deltaMasses;
				allIntensities[i]=intensities;
				retentionTimes[i]=scan.getScanStartTime();
			}

			int movingAverageLength=8; // expected points across the peak
			float[][] chromatograms=General.transposeMatrix(allIntensities);
			float[][] deltaMassByIon=General.transposeMatrix(allDeltaMasses);
			ArrayList<float[]> chromatogramList=new ArrayList<float[]>();
			ArrayList<FragmentIon> foundIons=new ArrayList<>();
			ArrayList<float[]> deltaMassList=new ArrayList<float[]>();
			for (int j = 0; j < chromatograms.length; j++) {
				if (primaryIonObjects[j].getIndex()>2&&General.sum(chromatograms[j])>0.0f) {
					if (sgSmoothBox.isSelected()) {
						chromatograms[j]=SkylineSGFilter.paddedSavitzkyGolaySmooth(chromatograms[j]);
					}
					if (backgroundSubtractBox.isSelected()) {
						chromatograms[j]=BackgroundSubtractionFilter.backgroundSubtractMovingMedian(chromatograms[j], movingAverageLength*10);
					}
					chromatogramList.add(chromatograms[j]);
					foundIons.add(primaryIonObjects[j]);
					deltaMassList.add(deltaMassByIon[j]);
				}
			}
			primaryIonObjects=foundIons.toArray(new FragmentIon[0]);
			
			TransitionRefinementData data;
			if (rtRange!=null&&rtRange.getRange()>0.0f) {
				data=TransitionRefiner.identifyTransitionsFromRTRange(entry.getPeptideModSeq(), entry.getPrecursorCharge(), entry.getRetentionTimeInSec(), 
						primaryIonObjects, chromatogramList, retentionTimes, rtRange, parameters);
			} else {
				data=TransitionRefiner.identifyTransitions(entry.getPeptideModSeq(), entry.getPrecursorCharge(), entry.getRetentionTimeInSec(), 
						primaryIonObjects, chromatogramList, retentionTimes, false, parameters);
			}
			fragmentTraces=getTraces(primaryIonObjects, data.getChromatograms(), data.getCorrelationArray(), retentionTimes, data.getRange(), -Float.MAX_VALUE);
			
			float minCorrelation=getMinimumCorrelation(data.getCorrelationArray());
			ArrayList<CategoricalData> deltaMassByIonList=new ArrayList<CategoricalData>();
			for (int j = 0; j < primaryIonObjects.length; j++) {
				if (data.getCorrelationArray()[j]>=minCorrelation) {
					TFloatArrayList deltaMasses=new TFloatArrayList();
					float[] deltaMassArray=deltaMassList.get(j);
					for (int i = 0; i < deltaMassArray.length; i++) {
						if (!Float.isNaN(deltaMassArray[i])&&data.getRange().contains(retentionTimes[i])) {
							deltaMasses.add(deltaMassArray[i]);
						}
					}
					deltaMassByIonList.add(new CategoricalData(foundIons.get(j).toString(), deltaMasses.toArray(), foundIons.get(j).getColor()));
				}				
			}
			
			final ChartPanel chartPanel = getChromatogramChartPanel(entry, fragmentTraces, precursorTraces, data.getRange());
	        
			dataPanel.add(chartPanel);
			
			JPanel rightInfoPanel=new JPanel(new GridLayout(2, 0));
			rightInfoPanel.setBackground(Color.WHITE);
			dataPanel.add(rightInfoPanel);
			
			MassTolerance fragmentTolerance = parameters.getFragmentTolerance();
			String deltaMassAxis=fragmentTolerance.isRelativeTolerance()?"Delta Mass (PPM)":"Delta Mass (AMU)";
			
			ChartPanel deltaMassPanel=Charter.getBoxplotChart("Delta Mass", "Ions", deltaMassAxis, 16, 16, deltaMassByIonList.toArray(new CategoricalData[0]), true);
			ValueAxis axis=deltaMassPanel.getChart().getCategoryPlot().getRangeAxis();
			axis.setRange(-fragmentTolerance.getToleranceThreshold(), fragmentTolerance.getToleranceThreshold());
			rightInfoPanel.add(deltaMassPanel);
			
			if (reference!=null) {
				ArrayList<LibraryEntry> references=reference.getEntries(entry.getPeptideModSeq(), entry.getPrecursorCharge(), false);
				if (references.size()==0) {
					references=reference.getEntries(entry.getLegacyPeptideModSeq(), entry.getPrecursorCharge(), false);
				}
				
				if (references.size()>0) {
					LibraryEntry ref=references.get(0);
					LibraryEntry acq=data.getEntry(ref, parameters);
	
					LibraryEntry butterfly=FragmentIonConsistencyCharter.getButterfly(acq, ref);
					ChartPanel chartPanelButterfly = Charter.getChart(new AnnotatedLibraryEntry(butterfly, parameters, true));

					Font font=new Font(Charter.BASE_FONT_NAME, Font.PLAIN, 18);
					XYTextAnnotation acquiredAnnotation = new XYTextAnnotation("Acquired", 10.0, 1.0);
					XYTextAnnotation libraryAnnotation = new XYTextAnnotation("Library", 10.0, -1.0);
					acquiredAnnotation.setTextAnchor(TextAnchor.TOP_LEFT);
					libraryAnnotation.setTextAnchor(TextAnchor.CENTER_LEFT);
					acquiredAnnotation.setPaint(Color.black);
					acquiredAnnotation.setFont(font);
					libraryAnnotation.setPaint(Color.black);
					libraryAnnotation.setFont(font);
					chartPanelButterfly.getChart().getXYPlot().addAnnotation(acquiredAnnotation);
					chartPanelButterfly.getChart().getXYPlot().addAnnotation(libraryAnnotation);

					rightInfoPanel.add(chartPanelButterfly);
				}
			}

		} catch (DataFormatException sqle) {
			Logger.errorLine("Data Format Error reading raw files!");
			Logger.errorException(sqle);
		} catch (SQLException sqle) {
			Logger.errorLine("SQL Error reading raw files!");
			Logger.errorException(sqle);
		} catch (IOException ioe) {
			Logger.errorLine("IO Error reading raw files!");
			Logger.errorException(ioe);
		} catch (Exception e) {
			Logger.errorLine("General Error reading raw files!");
			Logger.errorException(e);
		}
		
		mainSplit.setRightComponent(dataPanel);
		
		mainSplit.setDividerLocation(location);
	}

	private ChartPanel getChromatogramChartPanel(final InteractivePeptidePrecursor entry,
			ArrayList<XYTrace> fragmentTraces, ArrayList<XYTrace> precursorTraces, Range rtRange) {
		double globalMaxYFragment=0.0;
		double globalMaxYPrecursor=0.0;
		for (XYTrace xyTrace : precursorTraces) {
			if (xyTrace.getType()==GraphType.line) {
				globalMaxYPrecursor=Math.max(globalMaxYPrecursor, xyTrace.getMaxY());
			}
		}
		for (XYTrace xyTrace : fragmentTraces) {
			if (xyTrace.getType()==GraphType.boldline||xyTrace.getType()==GraphType.bolddashedline) {
				globalMaxYFragment=Math.max(globalMaxYFragment, xyTrace.getMaxY());
			}
		}
		
		fragmentTraces.add(new XYTrace(new float[] {rtRange.getStart()/60f, rtRange.getStop()/60f}, new float[] {(float)globalMaxYFragment, (float)globalMaxYFragment}, GraphType.area, "Boundaries", new Color(102, 204, 255, 50), 4.0f));
		precursorTraces.add(new XYTrace(new float[] {rtRange.getStart()/60f, rtRange.getStop()/60f}, new float[] {(float)globalMaxYPrecursor, (float)globalMaxYPrecursor}, GraphType.area, "Boundaries", new Color(102, 204, 255, 50), 4.0f));
		
		final ChartPanel chartPanel=ChromatogramCharter.createChart(Optional.ofNullable(precursorTraces), Optional.ofNullable(fragmentTraces), globalMaxYPrecursor, globalMaxYFragment);
		chartPanel.setMouseZoomable(false, false);

		CrosshairOverlay crosshairOverlay = new CrosshairOverlay();
		final Crosshair xCrosshair=new Crosshair(Double.NaN, Color.GRAY, new BasicStroke(0f));
		final EditableXYPoint zoomPoint=new EditableXYPoint();
		xCrosshair.setLabelVisible(true);
		crosshairOverlay.addDomainCrosshair(xCrosshair);
		chartPanel.addOverlay(crosshairOverlay);
		chartPanel.mouseDragged(null);

		chartPanel.getChart().setTitle(entry.getPeptideSeq());
		
		chartPanel.addMouseMotionListener(new MouseMotionListener() {
			@Override
			public void mouseMoved(MouseEvent e) {
		        Rectangle2D area=chartPanel.getScreenDataArea(e.getX(), e.getY());
		        if (area!=null) {
		        XYPlot plot = (XYPlot) chartPanel.getChart().getPlot();
		        ValueAxis xAxis = plot.getDomainAxis();
		        	double x = xAxis.java2DToValue(e.getX(), area, RectangleEdge.BOTTOM);
		        	xCrosshair.setValue(x);
		        }
			}
			
			@Override
			public void mouseDragged(MouseEvent e) {
		        Rectangle2D area=chartPanel.getScreenDataArea(e.getX(), e.getY());
		        if (area!=null) {
			        XYPlot plot = (XYPlot) chartPanel.getChart().getPlot();
			        ValueAxis xAxis = plot.getDomainAxis();
			        double x = xAxis.java2DToValue(e.getX(), area, RectangleEdge.BOTTOM);
			        xCrosshair.setValue(x);
			        
			        double prevX=zoomPoint.getX();
			        double prevY=zoomPoint.getY();
			        if (!Double.isNaN(prevX)&&!Double.isNaN(prevY)) {
				        Line2D zoomLine=new Line2D.Double(prevX, prevY, e.getX(), prevY);

				        Graphics2D g2 = (Graphics2D) chartPanel.getGraphics();
				        g2.setPaint(Color.gray);
				        g2.draw(zoomLine);
			        }
		        }
			}
		});
		
		chartPanel.addMouseListener(new MouseListener() {
			@Override
			public void mouseReleased(MouseEvent e) {
				double prevX=zoomPoint.getX();
				if (!Double.isNaN(prevX)) {
		            Rectangle2D area=chartPanel.getScreenDataArea(e.getX(), e.getY());
		        	double x=Math.max(area.getMinX(), Math.min(e.getX(), area.getMaxX()));
		        	
		        	double first=Math.min(x, prevX);
		        	double second=Math.max(x, prevX);
		        	
			        XYPlot plot = (XYPlot) chartPanel.getChart().getPlot();
			        ValueAxis xAxis = plot.getDomainAxis();
			        double plotX1 = xAxis.java2DToValue(first, area, RectangleEdge.BOTTOM);
			        double plotX2 = xAxis.java2DToValue(second, area, RectangleEdge.BOTTOM);
			        
			        entry.setRtRangeInSecs(new Range(plotX1*60.0, plotX2*60.0));
			        resetPeptide(entry);
		        }
			}
			
			@Override
			public void mousePressed(MouseEvent e) {
		        Rectangle2D area=chartPanel.getScreenDataArea(e.getX(), e.getY());
		        if (area!=null) {
		        	zoomPoint.setX(Math.max(area.getMinX(), Math.min(e.getX(), area.getMaxX())));
		        	zoomPoint.setY(Math.max(area.getMinY(), Math.min(e.getY(), area.getMaxY())));
		        }
		        else {
		        	zoomPoint.setX(null);
		        	zoomPoint.setY(null);
		        }
			}
			
			@Override
			public void mouseExited(MouseEvent e) {
			}
			
			@Override
			public void mouseEntered(MouseEvent e) {
			}
			
			@Override
			public void mouseClicked(MouseEvent e) {
			}
		});
		return chartPanel;
	}

	public Triplet<double[], float[], float[]> extract(Spectrum spectrum, FragmentIon[] ions, SearchParameters parameters) {
		MassTolerance acquiredTolerance=parameters.getFragmentTolerance();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();

		double[] actualTargetMasses=new double[ions.length];
		float[] actualTargetIntensities=new float[ions.length];
		float[] actualDeltaMasses=new float[ions.length];
		Arrays.fill(actualDeltaMasses, Float.NaN);
		
		for (int i = 0; i < ions.length; i++) {
			FragmentIon target=ions[i];
		
			int[] indicies=acquiredTolerance.getIndicies(acquiredMasses, target.getMass());
			float intensity=0.0f;
			float bestPeakIntensity=-1.0f;
			double bestPeakMass=0.0;
			
			for (int j=0; j<indicies.length; j++) {
				intensity+=acquiredIntensities[indicies[j]];
				
				if (acquiredIntensities[indicies[j]]>bestPeakIntensity) {
					bestPeakIntensity=acquiredIntensities[indicies[j]];
					bestPeakMass=acquiredMasses[indicies[j]];
				}
			}
			actualTargetIntensities[i]=intensity;
			actualTargetMasses[i]=bestPeakMass;
			if (intensity>0.0f&&bestPeakMass>0.0) {
				if (acquiredTolerance.isRelativeTolerance()) {
					// PPM
					actualDeltaMasses[i]=(float)(1000000*(bestPeakMass-target.getMass())/target.getMass());
				} else {
					actualDeltaMasses[i]=(float)(bestPeakMass-target.getMass());
				}
			}
		}
		
		return new Triplet<double[], float[], float[]>(actualTargetMasses, actualDeltaMasses, actualTargetIntensities);
	}

	private static final int minNumIons=6;
	private static ArrayList<XYTrace> getTraces(FragmentIon[] ions, ArrayList<float[]> chromatograms, float[] correlationArray, float[] rts, Range rtRange, float correlationThreshold) {
		ArrayList<XYTrace> xytraces=new ArrayList<XYTrace>();
		float minCorrelation = getMinimumCorrelation(correlationArray);
		
		for (int i=0; i<chromatograms.size(); i++) {
			if (correlationArray[i]<correlationThreshold) continue;
			
			float[] fs=chromatograms.get(i);
			
			TFloatArrayList rtSelectedIntensities=new TFloatArrayList();
			TFloatArrayList rtSelectedRTs=new TFloatArrayList();
			for (int j = 0; j < rts.length; j++) {
				if (rtRange.contains(rts[j])) {
					rtSelectedRTs.add(rts[j]);
					rtSelectedIntensities.add(fs[j]);
				}
			}
			
			Color c=ions[i].getColor();
			GraphType graphtype;
			GraphType backgroundgraphtype;
			float thickness;
			if (correlationArray[i]>=minCorrelation) {
				graphtype=GraphType.boldline;
				backgroundgraphtype=GraphType.dashedline;
				thickness=3.0f;
//			} else if (correlationArray[i]>TransitionRefiner.identificationCorrelationThreshold) {
//				graphtype=GraphType.boldline;
//				backgroundgraphtype=GraphType.line;
//				thickness=3.0f;
			} else {
				c=new Color(128, 128, 128, 128);
				graphtype=GraphType.bolddashedline;
				backgroundgraphtype=GraphType.dashedline;
				thickness=2.0f;
			}

			xytraces.add(TransitionRefiner.toXYTrace(rtSelectedIntensities.toArray(), rtSelectedRTs.toArray(), ""+i, c, rtRange, graphtype, thickness));

			if (rtRange!=null) {
				xytraces.add(TransitionRefiner.toXYTrace(fs, rts, ""+i, new Color(128, 128, 128, 128), null, backgroundgraphtype, thickness));
			}
		}
		return xytraces;
	}

	private static float getMinimumCorrelation(float[] correlationArray) {
		if (correlationArray.length==0) return TransitionRefiner.quantitativeCorrelationThreshold;
		
		float[] correlationClone=correlationArray.clone();
		Arrays.sort(correlationClone);
		float minCorrelation=correlationClone[Math.max(0,correlationClone.length-minNumIons)];
		return Math.min(TransitionRefiner.quantitativeCorrelationThreshold, minCorrelation);
	}

}
