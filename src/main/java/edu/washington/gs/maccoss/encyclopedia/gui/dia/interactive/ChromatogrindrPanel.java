package edu.washington.gs.maccoss.encyclopedia.gui.dia.interactive;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FileDialog;
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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;
import java.util.zip.DataFormatException;

import javax.swing.BorderFactory;
import javax.swing.Box;
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
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.RowFilter;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.panel.CrosshairOverlay;
import org.jfree.chart.plot.Crosshair;
import org.jfree.chart.plot.XYPlot;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.TextAnchor;

import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeptideXYPoint;
import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.NonstandardAminoAcidException;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.parameters.InstrumentSpecificSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaLibraryPredictionClient;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaPrecursor;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2023timsTOFModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.ims.AlphaPeptDeepIMSModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.ChronologerModel;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.FragmentIonConsistencyCharter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.PercentageLayout;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerProgress;
import edu.washington.gs.maccoss.encyclopedia.gui.massspec.ChromatogramCharter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.CategoricalData;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.EditableXYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BackgroundSubtractionFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

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
	
	private final JTable peptideTable;
	private final TableRowSorter<TableModel> rowSorter;
	private final PeptidePrecursorTableModel peptideModel;
	private final JTextField jtfFilter;
	private final JCheckBox jtfNotFilter=new JCheckBox("NOT");
	private final JComboBox<InstrumentSpecificSearchParameters> instrumentCombo=new JComboBox<InstrumentSpecificSearchParameters>(InstrumentSpecificSearchParameters.INSTRUMENTS);

	private final JCheckBox fragmentBox;
	private final JCheckBox sgSmoothBox;
	private final JCheckBox backgroundSubtractBox;
	private final SpinnerModel smallestIonNumbers=new SpinnerNumberModel(3, 0, 10, 1);

	private LibraryInterface reference=null;
	private StripeFileInterface dia=null;
	
	private final TObjectFloatHashMap<String> pastedRTs=new TObjectFloatHashMap<String>();
	private final TObjectFloatHashMap<String> libraryRTs=new TObjectFloatHashMap<String>();

	private final KoinaLibraryPredictionClient client;

	public static void main(String[] args) {
		File rawFile=new File("/Users/searleb/Documents/encyclopedia/small_file/bcs_2020jan16_hela_clib_3.mzML");
		File libraryFile=new File("/Users/searleb/Documents/encyclopedia/small_file/pan_human_library.dlib");
		final ChromatogrindrPanel browser=new ChromatogrindrPanel();
		launchBrowserPanel(browser);
		
		//browser.updateLibrary(libraryFile);
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
	}

	public static void launchBrowserPanel(final ChromatogrindrPanel browser) {
		final JFrame dialog=new JFrame("Chromatogrind Reanalysis Browser");

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
		
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(new Prosit2023timsTOFModel());
		models.add(new ChronologerModel());
		models.add(new AlphaPeptDeepIMSModel());
		this.client=new KoinaLibraryPredictionClient(models);
		
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

		fragmentBox=new JCheckBox("Fragments");
		fragmentBox.setSelected(true);
		fragmentBox.addActionListener(new ActionListener() {
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
		
		JSpinner ionSpinner = new JSpinner(smallestIonNumbers);
		ionSpinner.addChangeListener(new ChangeListener() {
			
			@Override
			public void stateChanged(ChangeEvent e) {
				updateToSelectedPeptide();
			}
		});
		LabeledComponent ionSpinnerLabel = new LabeledComponent("Min Ion", "The smallest ion index considered, 0 for unfiltered", ionSpinner);
		ionSpinnerLabel.setBackground(null);

		JPanel checkboxes=new JPanel(new FlowLayout());
		options.add(checkboxes);
		checkboxes.add(fragmentBox);
		checkboxes.add(Box.createHorizontalStrut(10));
		checkboxes.add(sgSmoothBox);
		checkboxes.add(Box.createHorizontalStrut(10));
		checkboxes.add(backgroundSubtractBox);
		checkboxes.add(Box.createHorizontalStrut(10));
		checkboxes.add(ionSpinnerLabel);
		
		JPanel buttons=new JPanel(new FlowLayout());
		options.add(buttons);
		
		JPanel in=new JPanel(new FlowLayout());
		int width=5;
		in.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(width, width, width, width), BorderFactory.createTitledBorder("In:")));
		JPanel out=new JPanel(new FlowLayout());
		out.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(width, width, width, width), BorderFactory.createTitledBorder("Out:")));
		buttons.add(in);
		buttons.add(Box.createHorizontalStrut(10));
		buttons.add(out);
		
		JButton importButton=new JButton("Import");
		importButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				importLibrary();
			}
		});
		in.add(importButton);
		
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
		in.add(pasteButton);
		
		JButton copyButton = new JButton("Copy");
		copyButton.addActionListener(new ActionListener() {	
			@Override
			public void actionPerformed(ActionEvent e) {
				copyTable();
			}
		});
		out.add(copyButton);
		
		JButton exportButton=new JButton("Export");
		exportButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				exportLibrary();
			}
		});
		out.add(exportButton);
		
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
	
	HashMap<String, AnnotatedLibraryEntry> entryMap=new HashMap<String, AnnotatedLibraryEntry>();

	public void pasteTable(String clip) {
		peptideModel.paste(clip);
		
		updatePredictions();
		
		pastedRTs.clear();
		for (InteractivePeptidePrecursor peptide : peptideModel.getAllEntries()) {
			pastedRTs.put(peptide.getPeptideModSeq(), peptide.getRetentionTimeInSec()/60f);
		}
		if (peptideTable.getRowCount()>0) {
			peptideTable.setRowSelectionInterval(0, 0);
		}
		peptideTable.requestFocus();
	}

	public void updatePredictions() {
		SearchParameters parameters=getParameters();
		ArrayList<KoinaPrecursor> pred=new ArrayList<KoinaPrecursor>();
		for (InteractivePeptidePrecursor entry : peptideModel.getAllEntries()) {
			if (!entryMap.containsKey(getPeptideIdentifier(entry))) {
				try {
					AminoAcidEncoding[] encoding = AminoAcidEncoding.getAAs(entry.getPeptideModSeq(), parameters.getAAConstants());
					boolean passes=client.checkPeptide(encoding, entry.getPrecursorCharge());
				
					KoinaPrecursor precursor = new KoinaPrecursor(encoding, 33f, entry.getPrecursorCharge());
					if (passes) {
						pred.add(precursor);
					}
				} catch (NonstandardAminoAcidException e) {
					Logger.errorLine("FAILED: "+entry.getPeptideModSeq());
				}
			}
		}
		int count=client.generatePredictions(KoinaLibraryPredictionClient.HTTPS_KOINA_WILHELMLAB_ORG_443, pred, new EmptyProgressIndicator(true));
		Logger.logLine("Generated "+count+" predictions");

		libraryRTs.clear();
		for (KoinaPrecursor precursor : pred) {
			AnnotatedLibraryEntry entry=precursor.toEntry(parameters.getAAConstants(), parameters);
			entryMap.put(getPeptideIdentifier(entry), entry);

			libraryRTs.put(entry.getPeptideModSeq(), entry.getRetentionTimeInSec()/60f);
		}
	}

	private static String getPeptideIdentifier(PeptidePrecursor entry) {
		String peptideModSeq = entry.getPeptideModSeq();
		byte precursorCharge = entry.getPrecursorCharge();
		return getPeptideIdentifier(peptideModSeq, precursorCharge);
	}

	private static String getPeptideIdentifier(String peptideModSeq, byte precursorCharge) {
		return peptideModSeq+":"+precursorCharge;
	}
	
	private Optional<LibraryEntry> getPrediction(PeptidePrecursor entry) {
		return Optional.ofNullable(entryMap.get(getPeptideIdentifier(entry)));
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
				
				libraryRTs.clear();
				for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants())) {
					libraryRTs.put(entry.getPeptideModSeq(), entry.getRetentionTimeInSec()/60f);
				}
				
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
	
	public void importLibrary() {
		FileDialog dialog=new FileDialog((JFrame)null, "Library files", FileDialog.LOAD);
		dialog.setFilenameFilter(new SimpleFilenameFilter(LibraryFile.DLIB, LibraryFile.ELIB, BlibFile.BLIB));
		dialog.setVisible(true);
		File[] fs=dialog.getFiles();

		if (fs==null&&fs.length==0&&fs[0]==null) {
			return;
		}
		File libraryFile=fs[0];

		try {
			LibraryFile library=new LibraryFile();
			library.openFile(libraryFile);
			ArrayList<LibraryEntry> entries=library.getAllEntries(false, new AminoAcidConstants());

			peptideModel.paste(entries);
			updatePredictions();
			
			pastedRTs.clear();
			for (InteractivePeptidePrecursor peptide : peptideModel.getAllEntries()) {
				pastedRTs.put(peptide.getPeptideModSeq(), peptide.getRetentionTimeInSec()/60f);
			}
			if (peptideTable.getRowCount()>0) {
				peptideTable.setRowSelectionInterval(0, 0);
			}
			peptideTable.requestFocus();

		} catch (DataFormatException dfe) {
			Logger.errorLine("Found Data Format error reading data from library file...");
			Logger.errorException(dfe);
		} catch (SQLException sqle) {
			Logger.errorLine("Found SQL error reading data from library file...");
			Logger.errorException(sqle);
		} catch (IOException ioe) {
			Logger.errorLine("Found IO error reading data from library file...");
			Logger.errorException(ioe);
		}
		
	}
	
	public void exportLibrary() {
		FileDialog dialog=new FileDialog((Frame)null, "Select a new or existing library file", FileDialog.SAVE);
		dialog.setFilenameFilter(new SimpleFilenameFilter(LibraryFile.DLIB));
		dialog.setVisible(true);
		File[] fs=dialog.getFiles();
		try {
			LibraryFile library=new LibraryFile();
			if (fs!=null&&fs.length>0&&fs[0]!=null) {
				File altFile=new File(fs[0].getParentFile(), fs[0].getName()+LibraryFile.DLIB);
				if (fs[0].exists()&&fs[0].canRead()) {
					Logger.logLine("Appending to existing file ["+fs[0].getName()+"].");
					library.openFile(fs[0]);
				} else if (altFile.exists()&&altFile.canRead()) {
						Logger.logLine("Appending to existing file ["+altFile.getName()+"].");
						library.openFile(altFile);
						fs[0]=altFile;
					
				} else {
					Logger.logLine("Creating new file ["+fs[0].getName()+"].");
					if (!fs[0].getName().toLowerCase().endsWith(LibraryFile.DLIB)) {
						File newFile=new File(fs[0].getAbsolutePath()+LibraryFile.DLIB);
						fs[0]=newFile;
					}
					library.openFile();
				}

				ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>(); 
						
				try {
					int lastUpdate=0;
					for (int i = 0; i < peptideModel.getRowCount(); i++) {
						int ratio=(100*i)/peptideModel.getRowCount();
						if (ratio>lastUpdate) {
							Logger.logLine(ratio+"% complete");
							lastUpdate=ratio;
						}
						InteractivePeptidePrecursor precursor=peptideModel.getSelectedRow(i);
						AnnotatedLibraryEntry entry=getLibraryEntry(precursor);
						if (entry!=null) {
							entries.add(entry);
						}
						
					}
				} catch (SQLException | IOException e) {
					Logger.errorException(e);
				}

				library.addEntries(entries, false);
				library.addProteinsFromEntries(entries);
				library.createIndices();

				library.saveAsFile(fs[0]);
				library.close();
				
				Logger.logLine("Finished exporting library!");
			}
		
		} catch (SQLException sqle) {
			Logger.errorLine("Found SQL error adding data to library file...");
			Logger.errorException(sqle);
		} catch (IOException ioe) {
			Logger.errorLine("Found IO error adding data to library file...");
			Logger.errorException(ioe);
		}
		
	}
	
	public AnnotatedLibraryEntry getLibraryEntry(final InteractivePeptidePrecursor entry) throws SQLException, IOException {
		SearchParameters parameters=getParameters();
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] primaryIonObjects=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false);
		if (((Integer)smallestIonNumbers.getValue())>1) {
			primaryIonObjects=FragmentIon.thresholdIons(primaryIonObjects, (Integer)smallestIonNumbers.getValue());
		}
		
		float rtInSec=entry.getRetentionTimeInSec();
		float minRTInSec = rtInSec-RT_EXTRACTION_MARGIN_IN_SEC;
		float maxRTInSec = rtInSec+RT_EXTRACTION_MARGIN_IN_SEC;

		Range rtRange=entry.getRTRange();
		if (rtRange!=null&&rtRange.getRange()>0.0f) {
			minRTInSec=rtRange.getStart();
			maxRTInSec=rtRange.getStop();
		}
		
		// get fragment traces
		ArrayList<FragmentScan> scans=dia.getStripes(entry.getPrecursorMZ(), minRTInSec, maxRTInSec, false);
		if (scans.size()==0) {
			// expand to full width
			scans=dia.getStripes(entry.getPrecursorMZ(), rtInSec-RT_EXTRACTION_MARGIN_IN_SEC, rtInSec+RT_EXTRACTION_MARGIN_IN_SEC, false);
		}
		if (scans.size()==0) {
			return null;
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
		
		AnnotatedLibraryEntry ref=FragmentationModel.generateEntry(entry.getPeptideModSeq(), dia.getFile().getName(), new HashSet<String>(), entry.getPrecursorCharge(), entry.getRetentionTimeInSec(), false, parameters);
		AnnotatedLibraryEntry acq=data.getEntry(ref, parameters);
		return acq;
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
		if (((Integer)smallestIonNumbers.getValue())>1) {
			primaryIonObjects=FragmentIon.thresholdIons(primaryIonObjects, (Integer)smallestIonNumbers.getValue());
		}
		
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
			ArrayList<PrecursorScan> trimmedPrecursors=new ArrayList<>();
			for (PrecursorScan spectrum : precursors) {
				if (entry.getPrecursorMZ()>spectrum.getIsolationWindowLower()&&entry.getPrecursorMZ()<spectrum.getIsolationWindowUpper()) {
					trimmedPrecursors.add(spectrum);
				}
			}
			Collections.sort(trimmedPrecursors);
			precursors=trimmedPrecursors;
			XYTraceInterface[] traceArray=ChromatogramExtractor.extractPrecursorChromatograms(parameters.getPrecursorTolerance(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), precursors, sgSmoothBox.isSelected(), backgroundSubtractBox.isSelected());
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
			
			
			final ChartPanel chartPanel = getChromatogramChartPanel(entry, fragmentTraces, precursorTraces, data.getRange());
	        
			dataPanel.add(chartPanel);
			
			JPanel rightInfoPanel=new JPanel(new GridLayout(0, 1));
			rightInfoPanel.setBackground(Color.WHITE);
			dataPanel.add(rightInfoPanel);

			// precursor delta mass
			PrecursorScanMap precursorMap=new PrecursorScanMap(precursors);
			MassTolerance precursorTolerance = parameters.getPrecursorTolerance();
			double precursorMz=parameters.getAAConstants().getChargedMass(entry.getPeptideModSeq(), entry.getPrecursorCharge());
			Pair<float[], TFloatArrayList[]> precursorIntegrations=precursorMap.integrateIsotopePacket(precursorMz, data.getRange(), entry.getPrecursorCharge(), parameters.getPrecursorTolerance());
			ArrayList<CategoricalData> deltaPrecursorMassByIonList=new ArrayList<CategoricalData>();
			for (int j = 0; j < precursorIntegrations.y.length; j++) {
				float[] deltaArray;
				if (precursorIntegrations.y[j].size()==0) {
					deltaArray=new float[] {0.0f};
				} else {
					deltaArray=precursorIntegrations.y[j].toArray();
				}
				deltaPrecursorMassByIonList.add(new CategoricalData(ChromatogramExtractor.getIsotopeShortName(j), deltaArray, ChromatogramExtractor.isotopeColors[j]));
			}
			
			String deltaPrecursorMassAxis=precursorTolerance.isRelativeTolerance()?"Delta Mass (PPM)":"Delta Mass (AMU)";

			PercentageLayout layout = new PercentageLayout(0, 2);
			layout.setColFraction(0, 0.33);
			layout.setColFraction(1, 0.67);
			JPanel deltaMassPanelGroup=new JPanel(layout);
			
			rightInfoPanel.add(deltaMassPanelGroup);
			
			ChartPanel deltaPrecursorMassPanel=Charter.getBoxplotChart("Delta Mass", "Ions", deltaPrecursorMassAxis, 16, 16, deltaPrecursorMassByIonList.toArray(new CategoricalData[0]), true);
			ValueAxis precursorAxis=deltaPrecursorMassPanel.getChart().getCategoryPlot().getRangeAxis();
			precursorAxis.setRange(-precursorTolerance.getToleranceThreshold(), precursorTolerance.getToleranceThreshold());
			deltaPrecursorMassPanel.getChart().setTitle("Precursor");
			deltaMassPanelGroup.add(deltaPrecursorMassPanel);
			
			// fragment delta mass
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
			
			MassTolerance fragmentTolerance = parameters.getFragmentTolerance();
			String deltaMassAxis=fragmentTolerance.isRelativeTolerance()?"Delta Mass (PPM)":"Delta Mass (AMU)";
			
			ChartPanel deltaMassPanel=Charter.getBoxplotChart("Delta Mass", "Ions", deltaMassAxis, 16, 16, deltaMassByIonList.toArray(new CategoricalData[0]), true);
			ValueAxis axis=deltaMassPanel.getChart().getCategoryPlot().getRangeAxis();
			axis.setRange(-fragmentTolerance.getToleranceThreshold(), fragmentTolerance.getToleranceThreshold());
			deltaMassPanel.getChart().setTitle("Fragments");
			deltaMassPanelGroup.add(deltaMassPanel);
			
			Optional<LibraryEntry> optionalRef=null;
			if (reference!=null) {
				ArrayList<LibraryEntry> references=reference.getEntries(entry.getPeptideModSeq(), entry.getPrecursorCharge(), false);
				if (references.size()==0) {
					references=reference.getEntries(entry.getLegacyPeptideModSeq(), entry.getPrecursorCharge(), false);
				}
				
				if (references.size()>0) {
					optionalRef=Optional.of(references.get(0));
				}
			}

			if (optionalRef==null) {
				optionalRef=getPrediction(entry);
			}
			
			if (optionalRef.isPresent()) {
				LibraryEntry ref=optionalRef.get();
				
				if ((Integer)smallestIonNumbers.getValue()>0) {
					ref=AnnotatedLibraryEntry.getAnnotationsOnly(ref, parameters, (Integer)smallestIonNumbers.getValue());
				}
				
				LibraryEntry acq=data.getEntry(ref, parameters);
				if ((Integer)smallestIonNumbers.getValue()>0) {
					acq=AnnotatedLibraryEntry.getAnnotationsOnly(acq, parameters, (Integer)smallestIonNumbers.getValue());
				}

				// THIS IS A SPEED ISSUE
				if (false&&pastedRTs.size()>0&&libraryRTs.size()>0) {
					Collection<XYPoint> points=new ArrayList<XYPoint>();
					pastedRTs.forEachEntry(new TObjectFloatProcedure<String>() {
						@Override
						public boolean execute(String a, float b) {
							if (libraryRTs.contains(a)) {
								float c=libraryRTs.get(a);
								points.add(new PeptideXYPoint(c, b, false, a));
							}
							return true;
						}
					});
					
					XYTrace backgroundRTs=new XYTrace(points, GraphType.bigpoint, "backgroundRTs", new Color(0f, 0f, 0f, 0.2f), 2.0f);
					XYTrace targetRT=new XYTrace(new float[] {ref.getScanStartTime()/60f}, new float[] {entry.getRetentionTimeInSec()/60f}, GraphType.bighollowpoint, "targetRTs", Color.red, 5.0f);
					
					ChartPanel rtPanel=Charter.getChart("Library Retention Time", "Acquired Retention Time", false, targetRT, backgroundRTs);
					rightInfoPanel.add(rtPanel);
				}

				PercentageLayout fraglayout = new PercentageLayout(0, 2);
				fraglayout.setColFraction(0, 0.33);
				fraglayout.setColFraction(1, 0.67);
				JPanel fragPanelGroup=new JPanel(fraglayout);
				rightInfoPanel.add(fragPanelGroup);

				// pad out for -1
				float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());
				predictedIsotopeDistribution=General.concatenate(new float[] {0f}, predictedIsotopeDistribution);
				
				ArrayList<XYTrace> precursorButterflyTraces=FragmentIonConsistencyCharter.getPrecursorButterfly(ChromatogramExtractor.isotopes, precursorIntegrations.x, predictedIsotopeDistribution, ChromatogramExtractor.isotopeColors);
				precursorButterflyTraces.add(new XYTrace(new float[] {-1.5f,  2.5f}, new float[] {0.0f, 0.0f}, GraphType.line, "base", Color.DARK_GRAY, 1f));
				ChartPanel chartPanelPrecursorButterfly = Charter.getChart("Precursor", "Relative Intensity", false, precursorButterflyTraces.toArray(new XYTrace[0]));
				((NumberAxis)chartPanelPrecursorButterfly.getChart().getXYPlot().getRangeAxis()).getAutoRangeIncludesZero();
				chartPanelPrecursorButterfly.getChart().getXYPlot().getDomainAxis().setRange(-1.5, 2.5);
				chartPanelPrecursorButterfly.getChart().getXYPlot().getRangeAxis().setRange(-1.15, 1.15);

				fragPanelGroup.add(chartPanelPrecursorButterfly);

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
				chartPanelButterfly.getChart().getXYPlot().getRangeAxis().setRange(-1.15, 1.15);
				chartPanelButterfly.getChart().getXYPlot().getRangeAxis().setLabel("Relative Intensity");
				chartPanelButterfly.getChart().setTitle((String)null);
				

				fragPanelGroup.add(chartPanelButterfly);
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
		
		if (!fragmentBox.isSelected()) {
			fragmentTraces=null;
		}
		
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
		            if (area==null) return;
		            
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
