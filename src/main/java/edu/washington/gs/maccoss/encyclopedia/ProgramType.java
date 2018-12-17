package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;

public enum ProgramType {
	EncyclopeDIA("EncyclopeDIA"), PecanPie("Walnut"), XCorDIA("XCorDIA"), CASiL("Thesaurus"), Global("Full EncyclopeDIA");
	
	private final String name;
	private ProgramType(String name) {
		this.name=name;
	}
	
	@Override
	public String toString() {
		return name;
	}
	
	public Version getVersion() {
		String version = ProgramType.class.getPackage().getImplementationVersion();
		Version localVersion=new Version(version);
		return localVersion;
	}
}
