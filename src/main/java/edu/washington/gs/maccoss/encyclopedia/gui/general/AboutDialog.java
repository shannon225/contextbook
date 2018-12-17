package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Desktop;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.net.URISyntaxException;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class AboutDialog {
	private static ImageIcon citationIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/citation.png"));
	public static void showAbout(JFrame parent, String shortName, String text, String citation, ImageIcon image) {
		final JDialog dialog=new JDialog(parent, "About "+shortName, true);
		
		JPanel p=new JPanel(new BorderLayout());
		
		p.add(new JLabel(image), BorderLayout.NORTH);
		p.add(new JEditorPane("text/html", "<html><center><p style=\"font-size:12px; font-family: Helvetica, sans-serif\">"+text), BorderLayout.CENTER);
		JPanel cite=new JPanel(new BorderLayout());
		JEditorPane citeTextBox=new JEditorPane("text/html", "<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\">"+citation);
		citeTextBox.setEditable(false);
		citeTextBox.addHyperlinkListener(new HyperlinkListener() {
		    public void hyperlinkUpdate(HyperlinkEvent e) {
		        if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
		            if (Desktop.isDesktopSupported()) {
		                try {
		                    Desktop.getDesktop().browse(e.getURL().toURI());
		                } catch (IOException e1) {
		                    Logger.errorException(e1);
		                } catch (URISyntaxException e1) {
		                    Logger.errorException(e1);
		                }
		            }
		        }
		    }
		});
		
		cite.add(citeTextBox, BorderLayout.CENTER);
		JLabel citeIcon=new JLabel(citationIcon);
		citeIcon.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		citeIcon.setBackground(Color.WHITE);
		citeIcon.setOpaque(true);
		cite.add(citeIcon, BorderLayout.WEST);
		p.add(cite, BorderLayout.SOUTH);
		p.setBackground(Color.WHITE);

		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				dialog.setVisible(false);
				dialog.dispose();
			}
		});

		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		buttons.add(okButton);
		buttons.setBackground(Color.WHITE);

		JPanel mainPane=new JPanel(new BorderLayout());
		mainPane.add(p, BorderLayout.CENTER);
		mainPane.add(buttons, BorderLayout.SOUTH);
		mainPane.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10), BorderFactory.createTitledBorder(shortName)));
		mainPane.setBackground(Color.WHITE);
		
		dialog.getContentPane().add(mainPane, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(400, 500);
		dialog.setVisible(true);
	}
}
