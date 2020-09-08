package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public class IndexedIonType {
	public final byte index;
	public final IonType type;
	public IndexedIonType(byte index, IonType type) {
		this.index=index;
		this.type=type;
	}
	public IndexedIonType(FragmentIon ion) {
		this.index=ion.getIndex();
		this.type=ion.getType();
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj==null||!(obj instanceof IndexedIonType)) return false;
		IndexedIonType o=(IndexedIonType)obj;
		return type==o.type&&index==o.index;
	}
	
	@Override
	public int hashCode() {
		// index * a prime larger than the largest likely index
		return index*101+type.hashCode();
	}	
}
