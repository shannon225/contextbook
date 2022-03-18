package edu.washington.gs.maccoss.encyclopedia.utils;

import junit.framework.TestCase;

import java.io.*;
import java.util.List;
import java.util.stream.Collectors;

public class LoggerTest extends TestCase {
    Throwable cause() {
        return new RuntimeException("Cause.");
    }

    Throwable exception() {
        return new RuntimeException("Exception.", cause());
    }

    /**
     * No assertions; check the console output!
     * We expect to see both the exception stacktrace AND it's cause.
     */
    public void testLogException() {
        Logger.logException(exception());
    }

    public void testErrorException() {
        Logger.errorException(exception());
    }

    public void testWriteStacktraceLines() throws IOException {
        try (ByteArrayOutputStream bos = new ByteArrayOutputStream()) {
            try (PrintStream ps = new PrintStream(bos)) {
                Logger.writeStacktraceLines(exception(), ps);
            }

            try (BufferedReader br = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(bos.toByteArray())))) {
                final List<String> lines =  br.lines().collect(Collectors.toList());

                assert lines.get(0).startsWith("[");
                assert lines.stream().anyMatch(l -> l.startsWith("Caused by:"));
                assert lines.stream().anyMatch(l -> l.contains("exception("));
                assert lines.stream().anyMatch(l -> l.contains("cause("));
            }
        }
    }
}