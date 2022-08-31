package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dialog;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DropTarget;
import java.awt.dnd.DropTargetDragEvent;
import java.awt.dnd.DropTargetDropEvent;
import java.awt.dnd.DropTargetEvent;
import java.awt.dnd.DropTargetListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EventObject;
import java.util.List;

import javax.swing.AbstractButton;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.event.CellEditorListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;

import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class FileChooserList extends JPanel {
	private static final long serialVersionUID=1L;
	//private static final String fileHeader = "File";
	private static final String deleteColumnHeader = " ";
	private static final ImageIcon fileDeleteIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/filedelete.png"));
	

	private final String fileType;
	private final FilenameFilter filter;
	
	private final ArrayList<File> files=new ArrayList<>();
	private final JButton chooseFile=new JButton("Add Files or Drag Into Box");
	private final FileTableModel model;
	private final JTable table;

    public static void main(String[] args) {
    	final JFrame f=new JFrame("TITLE");
		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});

		f.getContentPane().add(new FileChooserList(), BorderLayout.CENTER);

		f.pack();
		f.setSize(new Dimension(400,400));
		f.setVisible(true);
    }
    public FileChooserList() {
    	this("Any file type", new SimpleFilenameFilter(""));
    }
    
    public FileChooserList(String fileTypeLocal, FilenameFilter filterLocal) {
    	super(new BorderLayout());
    	
		this.fileType=fileTypeLocal;
		this.filter=filterLocal;
		
		this.model=new FileTableModel(fileTypeLocal);
		this.table=new JTable(model);

        // Create the drag and drop listener
        MyDragDropListener myDragDropListener = new MyDragDropListener();
        table.getColumn(fileTypeLocal).setCellRenderer(new PathCellRenderer());
        table.getColumn(deleteColumnHeader).setCellRenderer(new FileRemoverCellRendererAndEditor(table));
        table.getColumn(deleteColumnHeader).setCellEditor(new FileRemoverCellRendererAndEditor(table));
        table.setRowHeight(25);
        TableColumnModel columnModel = table.getColumnModel();
        columnModel.getColumn(0).setPreferredWidth(400);
        columnModel.getColumn(1).setPreferredWidth(25);
        
        JScrollPane pane=new JScrollPane(table);

        // Connect the label with a drag and drop listener
        new DropTarget(pane, myDragDropListener);
        
		chooseFile.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				askForFiles();
			}
		});

        // Add the label to the content
        this.add(pane, BorderLayout.CENTER);
        this.add(chooseFile, BorderLayout.NORTH);
	}
    
    public ArrayList<File> getFiles() {
		return files;
	}

	public void askForFiles() {
		File lastFile=null;
		if (files.size()>0) lastFile=files.get(files.size()-1);
		
		Component root=SwingUtilities.getRoot(FileChooserList.this);
		File[] fs;
		if (root instanceof Frame) {
			fs=FileChooserPanel.getFiles(lastFile, fileType, filter, (Frame)root, true);
		} else {
			fs=FileChooserPanel.getFiles(lastFile, fileType, filter, (Dialog)root, true);
			
		}
		files.addAll(Arrays.asList(fs));
    	notifyOfUpdate();
	}
	
    
    class FileTableModel extends AbstractTableModel {
		private static final long serialVersionUID = 1L;
		private final String fileHeader;
		

		public FileTableModel(String fileHeader) {
			super();
			this.fileHeader = fileHeader;
		}

		@Override
		public int getRowCount() {
			//if (files.size()==0) return 1;
			return files.size();
		}

		@Override
		public int getColumnCount() {
			return 2;
		}
		
		@Override
		public String getColumnName(int column) {
			if (column==0) return fileHeader;
			return deleteColumnHeader;
		}

		@Override
		public Object getValueAt(int rowIndex, int columnIndex) {
			if (columnIndex==0) {
				return files.get(rowIndex);
			}
			return "";
		}
		
		@Override
		public boolean isCellEditable(int rowIndex, int columnIndex) {
			if (columnIndex==1) return true;
			return false;
		}
    }
    
    class PathCellRenderer extends DefaultTableCellRenderer {
		private static final long serialVersionUID = 1L;

		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus,
				int row, int column) {
            JLabel c = (JLabel)super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
            if (value instanceof File) {
            	File file = (File)value;
				c.setToolTipText(file.getAbsolutePath());
				c.setText(file.getName());
            }
            return c;
        }
    }
    
	class FileRemoverCellRendererAndEditor implements TableCellRenderer, TableCellEditor {
		private JButton btn;
		private int row;

		FileRemoverCellRendererAndEditor(JTable table) {
			btn = new JButton(fileDeleteIcon);
			btn.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					File file=files.remove(row);
					notifyOfUpdate();

                	Logger.logLine("User deleted path is '" + file.getPath() + "'.");
				}
			});
		}

		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus,
				int row, int column) {
			return btn;
		}

		@Override
		public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row,
				int column) {
			this.row = row;
			return btn;
		}

		@Override
		public Object getCellEditorValue() {
			return true;
		}

		@Override
		public boolean isCellEditable(EventObject anEvent) {
			return true;
		}

		@Override
		public boolean shouldSelectCell(EventObject anEvent) {
			return true;
		}

		@Override
		public boolean stopCellEditing() {
			return true;
		}

		@Override
		public void cancelCellEditing() {
		}

		@Override
		public void addCellEditorListener(CellEditorListener l) {
		}

		@Override
		public void removeCellEditorListener(CellEditorListener l) {
		}
	}
    
    class MyDragDropListener implements DropTargetListener {

        @Override
        public void drop(DropTargetDropEvent event) {

            // Accept copy drops
            event.acceptDrop(DnDConstants.ACTION_COPY);

            // Get the transfer which can provide the dropped item data
            Transferable transferable = event.getTransferable();

            // Get the data formats of the dropped item
            DataFlavor[] flavors = transferable.getTransferDataFlavors();

            for (DataFlavor flavor : flavors) {
                try {
                    if (flavor.isFlavorJavaFileListType()) {
                        @SuppressWarnings("rawtypes")
						List list = (List) transferable.getTransferData(flavor);

                        for (Object obj : list) {
                        	if (obj instanceof File) {
	                        	File file=(File)obj;
	                        	if (filter.accept(file.getParentFile(), file.getName())) {
	                        		files.add(file);
	                        		notifyOfUpdate();
	                        		Logger.logLine("Added user dragged file path '" + file.getPath() + "'.");
	                        	} else {
		                            Logger.logLine("Could not user dragged add file path '" + file.getPath() + "'!");
	                        	}
                        	}
                        }

                    }

                } catch (Exception e) {
                    e.printStackTrace();

                }
            }

            event.dropComplete(true);

        }

        @Override
        public void dragEnter(DropTargetDragEvent event) {
        }

        @Override
        public void dragExit(DropTargetEvent event) {
        }

        @Override
        public void dragOver(DropTargetDragEvent event) {
        }

        @Override
        public void dropActionChanged(DropTargetDragEvent event) {
        }

    }

    public void addActionListener(ActionListener l) {
        listenerList.add(ActionListener.class, l);
    }

	private void notifyOfUpdate() {
		model.fireTableDataChanged();
		table.repaint();
		
		// Guaranteed to return a non-null array
        Object[] listeners = listenerList.getListenerList();
        ActionEvent e = null;
        // Process the listeners last to first, notifying
        // those that are interested in this event
        for (int i = listeners.length-2; i>=0; i-=2) {
            if (listeners[i]==ActionListener.class) {
                // Lazily create the event:
                if (e == null) {
					e = new ActionEvent(FileChooserList.this, ActionEvent.ACTION_PERFORMED, "Update file list");
                }
                ((ActionListener)listeners[i+1]).actionPerformed(e);
            }
        }
	}
}
