package edu.washington.gs.maccoss.encyclopedia.filereaders.mzml;


import com.google.common.base.Joiner;

import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class InstrumentComponentTranscoder {

	private enum Key {

		ORDER        ("order") {
			@Override
			public String encodePart(InstrumentComponent component) {
				return Integer.toString(component.order);
			}

			@Override
			protected InstrumentComponent.Builder decodeAppendImplementation(InstrumentComponent.Builder builder, String encodedEntryValue) {
				return builder.setOrder(Integer.parseInt(encodedEntryValue));
			}
		},
		CVREF        ("cvRef") {
			@Override
			public String encodePart(InstrumentComponent component) {
				return component.cvRef;
			}

			@Override
			protected InstrumentComponent.Builder decodeAppendImplementation(InstrumentComponent.Builder builder, String encodedEntryValue) {
				return builder.setCvRef(encodedEntryValue);
			}
		},
		ACCESSION_ID ("accessionId") {
			@Override
			public String encodePart(InstrumentComponent component) {
				return component.accessionId;
			}

			@Override
			protected InstrumentComponent.Builder decodeAppendImplementation(InstrumentComponent.Builder builder, String encodedEntryValue) {
				return builder.setAccessionId(encodedEntryValue);
			}
		},
		NAME         ("name") {
			@Override
			public String encodePart(InstrumentComponent component) {
				return component.name;
			}

			@Override
			protected InstrumentComponent.Builder decodeAppendImplementation(InstrumentComponent.Builder builder, String encodedEntryValue) {
				return builder.setName(encodedEntryValue);
			}
		},
		TYPE         ("type") {
			@Override
			public String encodePart(InstrumentComponent component) {
				return component.type.name;
			}

			@Override
			protected InstrumentComponent.Builder decodeAppendImplementation(InstrumentComponent.Builder builder, String encodedEntryValue) {
				InstrumentComponent.Type.getTypeByName(encodedEntryValue).ifPresent(builder::setType);
				return builder;
			}
		};

		final String key;

		Key(String key) {
			this.key = key;
		}

		public String encode(InstrumentComponent component) {
			return key + ":" + encodePart(component);
		}

		public abstract String encodePart(InstrumentComponent component);

		public InstrumentComponent.Builder decodeAppendGeneral(InstrumentComponent.Builder builder, String encodedEntries) {
			Pattern compile = Pattern.compile(key + ":([^" + ENTRY_DELIM + "]+)");
			Matcher matcher = compile.matcher(encodedEntries);
			if (matcher.find()) {
				decodeAppendImplementation(builder, matcher.group(1));
				return builder;
			} else {
				return builder;
			}
		}

		protected abstract InstrumentComponent.Builder decodeAppendImplementation(InstrumentComponent.Builder builder, String encodedEntryValue);

	}

	private static final String ENTRY_DELIM = ",";

	public static String encode(InstrumentComponent component) {
		return Joiner.on(ENTRY_DELIM).join(Arrays.stream(Key.values())
				.map(key -> key.encode(component))
				.iterator());
	}

	public static InstrumentComponent decode(String encoded) {
		InstrumentComponent.Builder builder = InstrumentComponent.builder();
		Arrays.stream(Key.values()).forEach(key -> {
			key.decodeAppendGeneral(builder, encoded);
		});
		return builder.build();
	}

}
