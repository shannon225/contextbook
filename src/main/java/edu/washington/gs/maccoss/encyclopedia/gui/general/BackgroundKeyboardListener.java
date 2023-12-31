package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.Component;
import java.awt.Container;
import java.awt.event.ContainerEvent;
import java.awt.event.ContainerListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.util.TreeSet;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.SearchGUIMain;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.interactive.ChromatogrindrPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class BackgroundKeyboardListener implements KeyListener, ContainerListener {
	TreeSet<EasterEgg> eggs=new TreeSet<EasterEgg>();
	
	public BackgroundKeyboardListener() {
		generateEasterEggs();
	}


	public void componentAdded(ContainerEvent e) {
		addKeyAndContainerListenerRecursively(e.getChild());
	}

	public void componentRemoved(ContainerEvent e) {
		removeKeyAndContainerListenerRecursively(e.getChild());
	}

	public void keyReleased(KeyEvent e) {
	}

	public void keyPressed(KeyEvent e) {
		int code=e.getKeyCode();

		if (e.isShiftDown()&&e.isControlDown()) {
			for (EasterEgg egg : eggs) {
				if (code==egg.getVirtualKeyCode()) {
					egg.run();
				}
			}
		}
	}

	public void keyTyped(KeyEvent e) {

	}

	public void addKeyAndContainerListenerRecursively(Component c) {
		c.addKeyListener(this);

		if (c instanceof Container) {
			Container cont=(Container)c;
			cont.addContainerListener(this);
			Component[] children=cont.getComponents();
			for (int i=0; i<children.length; i++) {
				addKeyAndContainerListenerRecursively(children[i]);
			}
		}
	}

	public void removeKeyAndContainerListenerRecursively(Component c) {
		c.removeKeyListener(this);

		if (c instanceof Container) {
			Container cont=(Container)c;
			cont.removeContainerListener(this);
			Component[] children=cont.getComponents();
			for (int i=0; i<children.length; i++) {
				removeKeyAndContainerListenerRecursively(children[i]);
			}
		}
	}

	void generateEasterEggs() {
		registerEasterEgg(new EasterEgg(KeyEvent.VK_H, "Show Easter Egg Help") {
			@Override void run() {
				Logger.logLine("Easter Egg Command Options:");
				for (EasterEgg egg : eggs) {
					Logger.logLine("\t"+egg.toString());
				}
			}
		});
		
		registerEasterEgg(new EasterEgg(KeyEvent.VK_A, "Launching with Advanced Options Enabled") {
			@Override void run() {
				Logger.logLine("Launching new window with Advanced Options Enabled");
				SearchGUIMain.runGUI(ProgramType.Global, true);
			}
		});
		
		registerEasterEgg(new EasterEgg(KeyEvent.VK_C, "Change Default Color") {
			@Override void run() {
				Logger.logLine("Changing default color");
				GUIParameters.requestUpdatedColor();
				Logger.logLine("Default color set to "+GUIParameters.getBaseColor().toString());
			}
		});
		
		registerEasterEgg(new EasterEgg(KeyEvent.VK_G, "Launch Interactive Chromatogram Visualizer") {
			@Override void run() {
				Logger.logLine("Launching Interactive Chromatogram Visualizer");
				ChromatogrindrPanel.launchBrowserPanel(new ChromatogrindrPanel());
			}
		});
	}

	void registerEasterEgg(EasterEgg egg) {
		eggs.add(egg);
	}

	abstract class EasterEgg implements Comparable<EasterEgg> {
		int virtualKeyCode;
		String name;

		public EasterEgg(int virtualKeyCode, String name) {
			this.virtualKeyCode=virtualKeyCode;
			this.name=name;
		}

		public int getVirtualKeyCode() {
			return virtualKeyCode;
		}

		@Override public String toString() {
			return KeyEvent.getKeyText(virtualKeyCode)+" --> "+name;
		}

		public int compareTo(EasterEgg o) {
			if (o==null) {
				return 1;
			} else {
				return (virtualKeyCode-o.virtualKeyCode);
			}
		}

		abstract void run();
	}

}
