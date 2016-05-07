package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.Color;

import javax.swing.JScrollPane;
import javax.swing.JTextPane;
import javax.swing.SwingUtilities;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.StyleConstants;

import edu.washington.gs.maccoss.encyclopedia.utils.LogRecorder;

public class LogConsole extends JScrollPane implements LogRecorder {
	private static final long serialVersionUID=1L;
	
	private final JTextPane text;
	private final Document doc;
	private final SimpleAttributeSet output;
	private final SimpleAttributeSet error;

	public LogConsole() {
		super();
		text=new JTextPane();
		text.setContentType("text/html");
		setViewportView(text);
		
		doc=text.getDocument();

		output=new SimpleAttributeSet();
		StyleConstants.setFontFamily(output, "SansSerif");
		StyleConstants.setFontSize(output, 10);

		error=new SimpleAttributeSet(output);
		StyleConstants.setBold(error, true);
		StyleConstants.setForeground(error, Color.red);
	}

	@Override
	public void log(String s) {
		log(s, output);
	}

	@Override
	public void logLine(String s) {
		log(s+"\n", output);
	}

	public void error(String s) {
		log(s, error);
	}

	@Override
	public void errorLine(String s) {
		log(s+"\n", error);
	}

	@Override
	public void logException(Throwable e) {
		log(throwableToString(e));
	}

	@Override
	public void errorException(Throwable e) {
		error(throwableToString(e));
	}

	public void log(final String s, final SimpleAttributeSet attributes) {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				try {
					doc.insertString(doc.getLength(), s, attributes);
				} catch (BadLocationException e) {
					throw new RuntimeException(e);
				}
				text.setCaretPosition(doc.getLength()-1);
			}
		});
	}

	public String throwableToString(Throwable e) {
		StringBuilder sb=new StringBuilder();
		sb.append(e.toString());
		sb.append("\n");
		for (StackTraceElement ste : e.getStackTrace()) {
			sb.append("\t");
			sb.append(ste.toString());
			sb.append("\n");
		}
		return sb.toString();
	}

}
