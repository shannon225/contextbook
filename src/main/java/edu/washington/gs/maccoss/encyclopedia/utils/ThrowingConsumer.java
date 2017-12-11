package edu.washington.gs.maccoss.encyclopedia.utils;

import java.util.function.Consumer;

@FunctionalInterface
public interface ThrowingConsumer<T> extends Consumer<T> {
	@Override
	default void accept(T t) {
		try {
			acceptThrows(t);
		} catch (final Exception e) {
			throw new RuntimeException(e);
		}
	}
	void acceptThrows(T t) throws Exception;
}
