package edu.washington.gs.maccoss.encyclopedia.utils;

public class EncyclopediaException extends RuntimeException {
	private static final long serialVersionUID=1L;

	public EncyclopediaException(String message, Throwable cause) {
		super(message, cause);
	}

	public EncyclopediaException(String message) {
		super(message);
	}

	public EncyclopediaException(Throwable cause) {
		super(cause);
	}

}
