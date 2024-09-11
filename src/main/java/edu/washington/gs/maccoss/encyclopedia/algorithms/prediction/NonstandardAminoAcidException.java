package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class NonstandardAminoAcidException extends EncyclopediaException {
	private static final long serialVersionUID = 1L;

	public NonstandardAminoAcidException(String message, Throwable cause) {
		super(message, cause);
	}

	public NonstandardAminoAcidException(String message) {
		super(message);
	}

	public NonstandardAminoAcidException(Throwable cause) {
		super(cause);
	}

}
