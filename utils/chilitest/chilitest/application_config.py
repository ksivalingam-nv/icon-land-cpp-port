# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

import argparse
import difflib
import inspect
import json
from pathlib import Path
from typing import Callable


class ConfigurationError(ValueError):
    pass


def validate_config(
    config: dict, schema: dict[str, tuple[Callable, set[str]]]
) -> None:
    """Validate an application config against the signatures it feeds.

    The config sections are splatted directly into the chilitest
    constructors (and test methods), so those signatures *are* the config
    schema; this check keeps it that way while turning the raw TypeError
    of a bad splat into messages that name the config section and key.

    `schema` maps a section path (e.g. "statistics", "test.z_score") to
    (target, injected): the callable the section is splatted into and the
    parameter names the application passes itself (not read from config).
    All problems are collected and reported in one go.
    """
    problems = []

    known_top_level = {path.split(".")[0] for path in schema}
    known_top_level.add("application_config_file")  # from the CLI overlay
    for name in sorted(config.keys() - known_top_level):
        problems.append(
            f"unknown section '{name}'" + _hint(name, known_top_level)
        )

    for path, (target, injected) in schema.items():
        section = config
        for part in path.split("."):
            section = section.get(part) if isinstance(section, dict) else None
            if section is None:
                problems.append(f"missing section '{path}'")
                break
        if section is None:
            continue
        if not isinstance(section, dict):
            problems.append(f"section '{path}' must be a json object")
            continue

        parameters = inspect.signature(target).parameters
        settable = {
            name
            for name, parameter in parameters.items()
            if name not in injected
            and name != "self"
            and parameter.kind
            not in (parameter.VAR_POSITIONAL, parameter.VAR_KEYWORD)
        }
        required = {
            name
            for name in settable
            if parameters[name].default is inspect.Parameter.empty
        }

        for name in sorted(section.keys() - settable):
            problems.append(
                f"section '{path}': unknown key '{name}'"
                + _hint(name, settable)
            )
        for name in sorted(required - section.keys()):
            problems.append(f"section '{path}': missing required key '{name}'")

    if problems:
        raise ConfigurationError(
            "Invalid application config:\n  - " + "\n  - ".join(problems)
        )


def _hint(name: str, candidates: set[str]) -> str:
    close_matches = difflib.get_close_matches(name, sorted(candidates), n=1)
    if close_matches:
        return f" (did you mean '{close_matches[0]}'?)"
    return f" (allowed: {sorted(candidates)})"


class ApplicationConfig:
    @classmethod
    def load(cls, application_config_file: str = None) -> dict:
        return cls._get_file_config(application_config_file)

    @classmethod
    def load_with_cli(cls, description: str = None) -> dict:
        """Load the application config named on the command line.

        Runtime information (paths, account, mkexp options) is given as
        regular command line options -- see `--help` -- and merged over the
        file config. The application-internal config keys are not exposed
        on the command line; everything else belongs in the config file.

        Note: not usable from a Jupyter notebook (the kernel occupies the
        command line); use `load` instead there.
        """
        cli_config, application_config_file = cls._get_cli_config(description)
        file_config = cls._get_file_config(application_config_file)

        # cli_config values override file_config values
        return cls._get_final_config(cli_config, file_config)

    @classmethod
    def _get_file_config(cls, application_config_file: str = None) -> dict:
        config = {}
        if application_config_file is not None:
            with open(application_config_file, "r") as f:
                config = json.load(f)
            config_dir = Path(application_config_file).resolve().parent
            for section_name, section in config.items():
                if isinstance(section, dict) and "include" in section:
                    config[section_name] = cls._resolve_include(
                        section, config_dir
                    )
        return config

    @classmethod
    def _resolve_include(cls, section: dict, config_dir: Path) -> dict:
        # A section may name a json file (path relative to the config file
        # declaring it) whose contents serve as defaults for the section.
        # This lets several application configs share one block -- e.g. the
        # analysis definition a test and its reference update must agree on
        # -- instead of duplicating it. Keys set in the section itself win.
        include_file = config_dir / section["include"]
        with open(include_file, "r") as f:
            included = json.load(f)
        if not isinstance(included, dict):
            raise ValueError(
                f"Included config must be a json object: {include_file}"
            )

        merged = {**included, **section}
        del merged["include"]
        return merged

    @classmethod
    def _get_cli_config(cls, description: str = None) -> tuple[dict, str]:
        # Only run-time *values* are settable on the command line; the
        # application-internal config keys they map onto are not exposed,
        # so they cannot be misspelled here. Analysis definitions and test
        # parameters belong in the application config file.
        parser = argparse.ArgumentParser(description=description)
        parser.add_argument(
            "application_config_file",
            help="application configuration json "
            "(analysis definition, test parameters, ...)",
        )
        parser.add_argument(
            "--mkexp-config",
            required=True,
            metavar="PATH",
            help="mkexp experiment configuration (<exp_name>.config)",
        )
        parser.add_argument(
            "--build-dir",
            required=True,
            metavar="DIR",
            help="ICON build directory (containing bin/icon)",
        )
        parser.add_argument(
            "--stats-dir",
            required=True,
            metavar="DIR",
            help="directory the reference statistics are written to/read from",
        )
        parser.add_argument(
            "--account",
            metavar="PROJECT",
            help="SLURM compute project "
            "(default: mkexp's fallback, e.g. ~/.acct)",
        )
        parser.add_argument(
            "--environment",
            metavar="NAME",
            help="mkexp environment, e.g. levante_gpu for a GPU build "
            "(default: chosen by mkexp for the host)",
        )
        parser.add_argument(
            "--output-dirname",
            metavar="NAME",
            help="name of an existing run_* output folder to analyse "
            "(default: found in the experiment work directory)",
        )
        parser.add_argument(
            "--mkexp-option",
            dest="mkexp_options",
            action="append",
            default=[],
            metavar="KEY=VALUE",
            help="additional mkexp option passed through to runscript "
            "generation (repeatable)",
        )
        args = parser.parse_args()

        file_info = {
            "mkexp_config_path": args.mkexp_config,
            "build_dir": args.build_dir,
            "stats_dir": args.stats_dir,
        }
        if args.output_dirname is not None:
            file_info["output_dirname"] = args.output_dirname

        options = {}
        for entry in args.mkexp_options:
            key, separator, value = entry.partition("=")
            if not separator or not key:
                parser.error(f"--mkexp-option expects KEY=VALUE, got '{entry}'")
            options[key] = value
        if args.account is not None:
            options["ACCOUNT"] = args.account
        if args.environment is not None:
            options["ENVIRONMENT"] = args.environment

        cli_config = {
            "file_info": file_info,
            "experiment": {"options": options},
        }
        return cli_config, args.application_config_file

    @classmethod
    def _get_final_config(cls, cli_config, file_config) -> dict:
        result = file_config.copy()

        for key, value in cli_config.items():
            if (
                key in result
                and isinstance(value, dict)
                and isinstance(result[key], dict)
            ):
                result[key] = cls._get_final_config(value, result[key])
            else:
                result[key] = value

        return result
