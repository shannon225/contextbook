package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public class AnnotatedFragmentIon extends FragmentIon {
	private final String annotation;

	public AnnotatedFragmentIon(double mass, String annotation) {
		super(mass, (byte)0, IonType.annotated);
		this.annotation = annotation;
	}
	
	@Override
	public String getName() {
		return annotation;
	}
}
