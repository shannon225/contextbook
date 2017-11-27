package edu.washington.gs.maccoss.encyclopedia.filereaders.mzml;

import com.google.common.base.Joiner;

import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class InstrumentIdTranscoder {

	private enum Key {
		CONFIGURATION_ID("configurationId") {
			@Override
			public String encodedPart(InstrumentId instrumentId) {
				return instrumentId.instrumentConfigurationId;
			}

			@Override
			protected InstrumentId.Builder decodeAppendImplementation(InstrumentId.Builder builder, String encodedEntryValue) {
				return builder.setInstrumentConfigurationId(encodedEntryValue);
			}
		},
		ACCESSION       ("accession") {
			@Override
			public String encodedPart(InstrumentId instrumentId) {
				return instrumentId.accession;
			}

			@Override
			protected InstrumentId.Builder decodeAppendImplementation(InstrumentId.Builder builder, String encodedEntryValue) {
				return builder.setAccession(encodedEntryValue);
			}
		},
		NAME            ("name") {
			@Override
			public String encodedPart(InstrumentId instrumentId) {
				return instrumentId.name;
			}

			@Override
			protected InstrumentId.Builder decodeAppendImplementation(InstrumentId.Builder builder, String encodedEntryValue) {
				return builder.setName(encodedEntryValue);
			}
		};

		private final String key;

		Key(String key) {
			this.key = key;
		}

		public abstract String encodedPart(InstrumentId instrumentId);

		public String encode(InstrumentId instrumentId) {
			return key + ":" + encodedPart(instrumentId);
		}

		public InstrumentId.Builder decodeAppendGeneral(InstrumentId.Builder builder, String encodedEntries) {
			Pattern compile = Pattern.compile(key + ":([^" + ENTRY_DELIM + "]+)");
			Matcher matcher = compile.matcher(encodedEntries);
			if (matcher.find()) {
				decodeAppendImplementation(builder, matcher.group(1));
			}
			return builder;
		}
		protected abstract InstrumentId.Builder decodeAppendImplementation(InstrumentId.Builder builder, String encodedEntryValue);
	}

	private static final String ENTRY_DELIM = ",";

	public static String encode(InstrumentId instrumentId) {
		return Joiner.on(ENTRY_DELIM).join(Arrays.stream(Key.values())
				.map(key -> key.encode(instrumentId))
				.iterator());

	}

	public static InstrumentId decode(String encoded) {
		InstrumentId.Builder builder = InstrumentId.builder();
		Arrays.stream(Key.values()).forEach(key -> {
			key.decodeAppendGeneral(builder, encoded);
		});
		return builder.build();
	}

}
